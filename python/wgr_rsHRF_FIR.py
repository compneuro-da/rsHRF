import numpy as np
from knee_pt import *
from scipy.sparse import lil_matrix
from scipy import stats

def wgr_BOLD_event_vector(N,matrix,thr,temporal_mask):
    """
    Detect BOLD event.
    event > thr & event < 3.1
    """
    data = lil_matrix((1,N))
    matrix = matrix[:,np.newaxis]
    if 0 in np.array(temporal_mask).shape:
        matrix = stats.zscore(matrix,ddof=1)
        for t in range(3,N-1):
            if ( matrix[t-1,0] > thr[0] and matrix[t-1,0] < thr[1] and np.all(matrix[t-3:t-1,0] < matrix[t-1,0]) and np.all(matrix[t-1,0] > matrix[t:t+2,0]) ):
                data[0,t-1] = 1
    else:
        datm = np.mean(matrix[temporal_mask])
        datstd = np.std(matrix[temporal_mask])
        datstd[datstd==0] = 1
        matrix = np.divide((matrix - datm),datstd)
        for t in range(3,N-1):
            if temporal_mask[t]:
                if ( matrix[t-1,0] > thr[0] and matrix[t-1,0] < thr[1] and np.all(matrix[t-3:t-1,0] < matrix[t-1,0]) and np.all(matrix[t-1,0] > matrix[t:t+2,0]) ):
                    data[0,t-1] = 1
    return data

def tor_make_deconv_mtx3(sf,tp,eres):
    docenter = 0
    if type(sf) is not dict:
        sf2 = {}
        for i in range(0,sf.shape[1]):
            sf2[i] = sf[:,i]
        sf = sf2
    if type(tp) is int:
        tp = np.tile(tp,(1,len(sf)))
    if len(tp) != len(sf):
        print('timepoints vectors (tp) and stick function (sf) lengths do not match!')
        return
    tbefore = 0
    nsess = len(sf)

    numtrs = int(np.around(np.amax(sf[0].shape) / eres))
    myzeros = np.zeros((numtrs,1))
    DX = np.zeros((numtrs,1))

    for i in range(0,len(sf)):
        Snumtrs = np.amax(sf[i].shape) / eres
        if(Snumtrs!=np.round(Snumtrs)):
            print('length not evenly divisible by eres')
        if(numtrs!=Snumtrs):
            print('different length than sf[0]')

        inums = np.nonzero(sf[i] > 0)[0]
        inums = inums / eres
        inums = np.ceil(inums).astype(int)
        sf[i] = np.ravel(myzeros)
        sf[i][inums] = 1

    index = 0
    for i in range(0,len(sf)):
        if tbefore!=0:
            for j in range(tbefore-1,-1,-1):
                sf_temp = sf[i][j:]
                sf_temp = sf_temp[:,np.newaxis]
                mysf = np.concatenate((sf_temp,np.zeros((j,1))))
                if index==0:
                    DX[:,index] = np.ravel(mysf)
                else:
                    DX = np.column_stack((DX,mysf))
                index += 1

        if index == 0:
            DX[:, index] = sf[i]
        else:
            DX = np.column_stack((DX,sf[i]))

        index += 1
        inums = np.nonzero(sf[i]==1)[0]

        for j in range(1,np.ravel(tp)[i]):
            myzeros = np.zeros((numtrs, 1))
            inums = inums + 1
            reg = myzeros
            inums = inums[inums<numtrs]
            reg[inums] = 1
            while (np.amax(reg.shape) < DX.shape[0]):
                reg = np.concatenate((reg,np.zeros(1,1)))
            DX = np.column_stack((DX, reg))
            index += 1

    if nsess < 2:
        DX = np.column_stack((DX,np.ones((DX.shape[0],1))))
    else:
        X = np.zeros((DX.shape[0],1))
        index = 0
        scanlen = DX.shape[0] / nsess
        if np.around(scanlen) != scanlen:
            print('Model length is not an even multiple of scan length.')
        for startimg in range(0,DX.shape[0],int(np.around(scanlen))):
            if index == 0:
                X[startimg:startimg+int(np.around(scanlen)),index] = 1
            else:
                X_temp = np.zeros((DX.shape[0],1))
                X_temp[startimg:startimg+int(np.around(scanlen)),0] = 1
                X = np.column_stack((X,X_temp))
            index += 1
        DX = np.column_stack((DX,X))

    if docenter:
        wh = np.arange(1,DX.shape[1]-nsess+1)
        DX[:,wh] = DX[:,wh] - np.tile(np.mean(DX[:,wh]),(DX.shape[0],1))
    return DX, sf

