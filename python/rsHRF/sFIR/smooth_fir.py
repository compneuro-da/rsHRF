import numpy as np
from scipy.sparse import lil_matrix
from scipy import stats
from joblib import Parallel, delayed
from joblib import load, dump
import tempfile
import shutil
import os
import warnings
from ..processing import knee

warnings.filterwarnings("ignore")

def wgr_BOLD_event_vector(N, matrix, thr, k, temporal_mask):
    """
    Detect BOLD event.
    event > thr & event < 3.1
    """
    data = lil_matrix((1, N))
    matrix = matrix[:, np.newaxis]
    if 0 in np.array(temporal_mask).shape:
        matrix = stats.zscore(matrix, ddof=1)
        matrix = np.nan_to_num(matrix)
        for t in range(1 + k, N - k + 1):
            if matrix[t - 1, 0] > thr[0] and \
                    np.all(matrix[t - k - 1:t - 1, 0] < matrix[t - 1, 0]) and \
                    np.all(matrix[t - 1, 0] > matrix[t:t + k, 0]):
                data[0, t - 1] = 1
    else:
        datm = np.mean(matrix[temporal_mask])
        datstd = np.std(matrix[temporal_mask])
        datstd[datstd == 0] = 1
        matrix = np.divide((matrix - datm), datstd)
        for t in range(1 + k, N - k + 1):
            if temporal_mask[t-1]:
                if matrix[t - 1, 0] > thr[0] and \
                        np.all(matrix[t - k - 1:t - 1, 0] < matrix[t - 1, 0]) and \
                        np.all(matrix[t - 1, 0] > matrix[t:t + k, 0]):
                    data[0, t - 1] = 1
    return data


def tor_make_deconv_mtx3(sf, tp, eres):
    docenter = 0
    if type(sf) is not dict:
        sf2 = {}
        for i in range(0, sf.shape[1]):
            sf2[i] = sf[:, i]
        sf = sf2
    if type(tp) is int:
        tp = np.tile(tp, (1, len(sf)))
    if len(tp) != len(sf):
        print('timepoints vectors (tp) and \
        stick function (sf) lengths do not match!')
        return
    tbefore = 0
    nsess = len(sf)

    numtrs = int(np.around(np.amax(sf[0].shape) / eres))
    myzeros = np.zeros((numtrs, 1))
    DX = np.zeros((numtrs, 1))

    for i in range(0, len(sf)):
        Snumtrs = np.amax(sf[i].shape) / eres
        if(Snumtrs != np.round(Snumtrs)):
            print('length not evenly divisible by eres')
        if(numtrs != Snumtrs):
            print('different length than sf[0]')

        inums = np.nonzero(sf[i] > 0)[0]
        inums = inums / eres
        inums = np.ceil(inums).astype(int)
        sf[i] = np.ravel(myzeros)
        sf[i][inums] = 1

    index = 0
    for i in range(0, len(sf)):
        if tbefore != 0:
            for j in range(tbefore - 1, -1, -1):
                sf_temp = sf[i][j:]
                sf_temp = sf_temp[:, np.newaxis]
                mysf = np.concatenate((sf_temp, np.zeros((j, 1))))
                if index == 0:
                    DX[:, index] = np.ravel(mysf)
                else:
                    DX = np.column_stack((DX, mysf))
                index += 1

        if index == 0:
            DX[:, index] = sf[i]
        else:
            DX = np.column_stack((DX, sf[i]))

        index += 1
        inums = np.nonzero(sf[i] == 1)[0]

        for j in range(1, np.ravel(tp)[i]):
            myzeros = np.zeros((numtrs, 1))
            inums = inums + 1
            reg = myzeros
            inums = inums[inums < numtrs]
            reg[inums] = 1
            while (np.amax(reg.shape) < DX.shape[0]):
                reg = np.concatenate((reg, np.zeros(1, 1)))
            DX = np.column_stack((DX, reg))
            index += 1

    if nsess < 2:
        DX = np.column_stack((DX, np.ones((DX.shape[0], 1))))
    else:
        X = np.zeros((DX.shape[0], 1))
        index = 0
        scanlen = DX.shape[0] / nsess
        if np.around(scanlen) != scanlen:
            print('Model length is not an even multiple of scan length.')
        for startimg in range(0, DX.shape[0], int(np.around(scanlen))):
            if index == 0:
                X[startimg:startimg + int(np.around(scanlen)), index] = 1
            else:
                X_temp = np.zeros((DX.shape[0], 1))
                X_temp[startimg:startimg + int(np.around(scanlen)), 0] = 1
                X = np.column_stack((X, X_temp))
            index += 1
        DX = np.column_stack((DX, X))

    if docenter:
        wh = np.arange(1, DX.shape[1] - nsess + 1)
        DX[:, wh] = DX[:, wh] - np.tile(np.mean(DX[:, wh]), (DX.shape[0], 1))
    return DX, sf


def Fit_sFIR2(tc, TR, Runs, T, mode):
    DX, sf = tor_make_deconv_mtx3(Runs, T, 1)
    DX2 = DX[:, 0:T]
    num = T

    if mode == 1:
        C = np.arange(1, num + 1).reshape((1, num)).conj().T\
            .dot(np.ones((1, num)))
        h = np.sqrt(1 / (7 / TR))

        v = 0.1
        sig = 1

        R = v * np.exp(-h / 2 * (C - C.conj().T) ** 2)
        RI = np.linalg.inv(R)

        b = np.linalg.solve((DX2.conj().T.dot(DX2) + sig ** 2 * RI),
                            DX2.conj().T).dot(tc)
        e = tc - DX2.dot(b)

    elif mode == 0:
        b = np.linalg.pinv(DX).dot(tc)
        e = tc - DX.dot(b)
        b = b[0:T]

    hrf = b
    return hrf, e


def wgr_rsHRF_FIR(data, para, temporal_mask):
    para['temporal_mask'] = temporal_mask
    N, nvar = data.shape
    if np.count_nonzero(para['thr']) == 1:
        para['thr'] = np.array([para['thr'], np.inf])

    folder = tempfile.mkdtemp()
    data_folder = os.path.join(folder, 'data')
    dump(data, data_folder)
    data = load(data_folder, mmap_mode='r')

    results = Parallel(n_jobs=-1)(delayed(wgr_FIR_estimation_HRF)(data, i, para, N) for i in range(0, nvar))
    beta_rshrf, event_bold = zip(*results)

    try:
        shutil.rmtree(folder)
    except:
        print("Failed to delete: " + folder)

    return np.array(beta_rshrf).T, np.array(event_bold)


def wgr_FIR_estimation_HRF(data, i, para, N):
    if para['estimation'] == 'sFIR':
        firmode = 1
    else:
        firmode = 0
    dat = data[:, i]

    if 'localK' not in para:
        if para['TR']<=2:
            localK = 1
        else:
            localK = 2
    else:
        localK = para['localK']

    u = wgr_BOLD_event_vector(N, dat, para['thr'], localK, para['temporal_mask'])
    u = u.toarray().flatten(1).ravel().nonzero()[0]

    lag = para['lag']
    nlag = np.amax(lag.shape)
    len_bin = int(np.floor(para['len'] / para['TR']))

    hrf = np.zeros((len_bin, nlag))
    Cov_E = np.zeros((1, nlag))
    kk = 0

    for i_lag in range(1, nlag + 1):
        RR = u - i_lag
        RR = RR[RR >= 0]
        if RR.size != 0:
            design = np.zeros((N, 1))
            design[RR] = 1
            hrf_kk, e3 = Fit_sFIR2(dat, para['TR'], design, len_bin, firmode)
            hrf[:, kk] = np.ravel(hrf_kk)
            Cov_E[:, kk] = np.cov(np.ravel(e3))
        else:
            Cov_E[:, kk] = np.inf
        kk += 1

    placeholder, ind = knee.knee_pt(np.ravel(Cov_E))
    rsH = hrf[:, ind + 1]
    return rsH, u
