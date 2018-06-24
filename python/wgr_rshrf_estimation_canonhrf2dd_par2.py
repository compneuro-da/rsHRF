import numpy as np
from spm import *
from scipy.sparse import lil_matrix
from scipy import stats


def wgr_rshrf_estimation_canonhrf2dd_par2(data, xBF, temporal_mask):
    N, nvar = data.shape
    bf = wgr_spm_get_canonhrf(xBF)
    bf2 = wgr_spm_Volterra(bf, xBF)
    if bf2 != []:
        bf = np.column_stack((bf, bf2))

    length = xBF['len']
    beta_hrf = []
    event_bold = []

    for i in range(nvar):
        beta_hrf_out, event_bold_out = \
            wgr_hrf_estimation_canon(data[:, i], xBF, length,
                                     N, bf, temporal_mask)
        beta_hrf.append(beta_hrf_out)
        event_bold.append(event_bold_out)
    return np.array(beta_hrf).T, bf, np.array(event_bold)


def wgr_spm_Volterra(bf, xBF):
    bf2 = []
    if 'Volterra' in xBF:
        if xBF['Volterra'] == 2:
            bf2 = []
            for p in range(bf.shape[1]):
                for q in range(bf.shape[1]):
                    bf2.append((bf[:, p] * bf[:, q]))
            bf2 = np.array(bf2).T
            bf2 = spm_orth(bf2)
    return bf2


def wgr_spm_get_canonhrf(xBF):
    dt = xBF['dt']
    fMRI_T = xBF['T']

    bf = spm_hrf(dt, P=None, fMRI_T=fMRI_T)
    p = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    p[len(p) - 1] = xBF['len']

    bf = spm_hrf(dt, p, fMRI_T)
    bf = bf[:, np.newaxis]

    if xBF['TD_DD']:
        dp = 1
        p[5] = p[5] + dp
        D = (bf[:, 0] - spm_hrf(dt, p, fMRI_T)) / dp
        D = D[:, np.newaxis]
        bf = np.append(bf, D, axis=1)
        p[5] = p[5] - dp
        if xBF['TD_DD'] == 2:
            dp = 0.01
            p[2] = p[2] + dp
            D = (bf[:, 0] - spm_hrf(dt, p, fMRI_T)) / dp
            D = D[:, np.newaxis]
            bf = np.append(bf, D, axis=1)

    bf = spm_orth(bf)
    return bf


def wgr_hrf_estimation_canon(dat, xBF, length, N, bf, temporal_mask):
    """
    Estimate HRF
    """
    thr = xBF['thr']
    u0 = wgr_BOLD_event_vector(N, dat, thr, temporal_mask)
    u = np.append(u0.toarray(), np.zeros((xBF['T'] - 1, N)), axis=0)
    u = np.reshape(u, (1, - 1), order='F')
    beta, lag = wgr_hrf_fit(dat, length, xBF, u, N, bf)
    beta_hrf = beta
    beta_hrf = np.append(beta_hrf, lag)
    return beta_hrf, u0.toarray()[0].nonzero()[0]


def wgr_BOLD_event_vector(N, matrix, thr, temporal_mask):
    """
    Detect BOLD event.
    event > thr & event < 3.1
    """
    data = lil_matrix((1, N))
    matrix = matrix[:, np.newaxis]
    if 0 in np.array(temporal_mask).shape:
        matrix = stats.zscore(matrix, ddof=1)
        matrix = np.nan_to_num(matrix)
        for t in range(3, N - 1):
            if thr < matrix[t - 1, 0] < 3.1 and \
                    np.all(matrix[t - 3:t - 1, 0] < matrix[t - 1, 0]) and \
                    np.all(matrix[t - 1, 0] > matrix[t:t + 2, 0]):
                data[0, t - 1] = 1
    else:
        datm = np.mean(matrix[temporal_mask])
        datstd = np.std(matrix[temporal_mask])
        datstd[datstd == 0] = 1
        matrix = np.divide((matrix - datm), datstd)
        for t in range(3, N - 1):
            if temporal_mask[t]:
                if thr < matrix[t - 1, 0] < 3.1 and \
                        np.all(matrix[t - 3:t - 1, 0] < matrix[t - 1, 0]) and \
                        np.all(matrix[t - 1, 0] > matrix[t:t + 2, 0]):
                    data[0, t - 1] = 1
    return data


def wgr_hrf_fit(dat, length, xBF, u, N, bf):
    """
    @u    - BOLD event vector (microtime).
    @nlag - time lag from neural event to BOLD event
    """
    lag = xBF['lag']
    AR_lag = xBF['AR_lag']
    nlag = len(lag)
    erm = np.zeros((1, nlag))
    beta = np.zeros((bf.shape[1] + 1, nlag))
    for i in range(nlag):
        u_lag = np.append(u[0, lag[i]:], np.zeros((1, lag[i]))).T
        erm[0, i], beta[:, i] = \
            wgr_glm_estimation(dat, u_lag, bf, xBF['T'], xBF['T0'], AR_lag)

    idx = np.argmin(erm)
    return beta[:, idx], lag[idx]


def wgr_glm_estimation(dat, u, bf, T, T0, AR_lag):
    """
    @u - BOLD event vector (microtime).
    """
    nscans = dat.shape[0]
    x = wgr_onset_design(u, bf, T, T0, nscans)
    X = np.append(x, np.ones((nscans, 1)), axis=1)
    res_sum, Beta = wgr_glsco(X, dat, AR_lag)
    return res_sum, Beta


def wgr_onset_design(u, bf, T, T0, nscans):
    """
    @u - BOLD event vector (microtime).
    @bf - basis set matrix
    @T - microtime resolution (number of time bins per scan)
    @T0 - microtime onset (reference time bin, see slice timing)
    """
    ind = np.arange(0, max(u.shape))
    X = np.empty((0, len(ind)))
    for p in range(bf.shape[1]):
        x = np.convolve(u, bf[:, p])
        x = x[ind]
        X = np.append(X, [x], axis=0)
    X = X.T
    """
    Resample regressors at acquisition times
    """
    X = X[(np.arange(0, nscans) * T) + (T0 - 1), :]
    return X


def wgr_glsco(X, Y, AR_lag=1, max_iter=20):
    """
    Linear regression when disturbance terms follow AR(p)
    -----------------------------------
    Model:
    Yt = Xt * Beta + ut ,
    ut = Phi1 * u(t-1) + ... + Phip * u(t-p) + et
    where et ~ N(0,s^2)
    -----------------------------------
    Algorithm:
    Cochrane-Orcutt iterated regression (Feasible generalized least squares)
    -----------------------------------
    Usage:
    Y = dependent variable (n * 1 vector)
    X = regressors (n * k matrix)
    AR_lag = number of lags in AR process
    -----------------------------------
    Returns:
    Beta = estimator corresponding to the k regressors
    """
    nobs, nvar = X.shape
    Beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    resid = Y - (X.dot(Beta))

    max_tol = min(1e-6, max(np.absolute(Beta)) / 1000)
    for r in range(max_iter):
        Beta_temp = Beta
        X_AR = np.zeros((nobs - (2 * AR_lag), AR_lag))

        for m in range(AR_lag):
            X_AR[:, m] = resid[AR_lag - m - 1:nobs - AR_lag - m - 1]

        Y_AR = resid[AR_lag:nobs - AR_lag]
        AR_para = np.linalg.lstsq(X_AR, Y_AR, rcond=None)[0]

        X_main = X[AR_lag:nobs, :]
        Y_main = Y[AR_lag:nobs]

        for m in range(AR_lag):
            X_main = \
                X_main - (AR_para[m] * (X[AR_lag - m - 1:nobs - m - 1, :]))
            Y_main = Y_main - (AR_para[m] * (Y[AR_lag - m - 1:nobs - m - 1]))

        Beta = np.linalg.lstsq(X_main, Y_main, rcond=None)[0]
        resid = Y[AR_lag:nobs] - X[AR_lag:nobs, :].dot(Beta)
        if(max(np.absolute(Beta - Beta_temp)) < max_tol):
            break

    res_sum = np.sum(resid ** 2)
    return res_sum, Beta
