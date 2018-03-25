import numpy as np
from spm import *
from scipy.sparse import lil_matrix
from scipy import stats

def wgr_onset_design(u,bf,T,T0,nscans):
    """
    @u - BOLD event vector (microtime).
    @bf - basis set matrix
    @T - microtime resolution (number of time bins per scan)
    @T0 - microtime onset (reference time bin, see slice timing)
    """
    ind = np.arange(0,max(u.shape))
    X = np.empty((0,len(ind)))
    for p in range(bf.shape[1]):
        x = np.convolve(u,bf[:,p])
        x = x[ind]
        X = np.append(X,[x],axis=0)
    X = X.T
    """
    Resample regressors at acquisition times
    """
    X = X[(np.arange(0,nscans)*T)+(T0-1),:]
    return X
