import numpy as np
from spm import *
from scipy.sparse import lil_matrix
from scipy import stats

def wgr_glm_estimation(dat,u,bf,T,T0,AR_lag):
    """
    @u - BOLD event vector (microtime).
    """
    nscans = dat.shape[0]
    x = wgr_onset_design(u,bf,T,T0,nscans)
    X = np.append(x,np.ones((nscans,1)),axis=1)
    res_sum, Beta = wgr_glsco(X,dat,AR_lag)
    return res_sum, Beta

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

def wgr_glsco(X,Y,AR_lag=1,max_iter=20):
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
    Beta = np.linalg.lstsq(X,Y)[0]
    resid = Y - (X.dot(Beta))

    max_tol = min(1e-6,max(np.absolute(Beta))/1000)
    for r in range(max_iter):
        Beta_temp = Beta
        X_AR = np.zeros((nobs-(2*AR_lag),AR_lag))

        for m in range(AR_lag):
            X_AR[:,m] = resid[AR_lag-m-1:nobs-AR_lag-m-1]

        Y_AR = resid[AR_lag:nobs-AR_lag]
        AR_para = np.linalg.lstsq(X_AR,Y_AR)[0]

        X_main = X[AR_lag:nobs,:]
        Y_main = Y[AR_lag:nobs]

        for m in range(AR_lag):
            X_main = X_main - (AR_para[m]*(X[AR_lag-m-1:nobs-m-1,:]))
            Y_main = Y_main - (AR_para[m]*(Y[AR_lag-m-1:nobs-m-1]))

        Beta = np.linalg.lstsq(X_main,Y_main)[0]
        resid = Y[AR_lag:nobs] - X[AR_lag:nobs,:].dot(Beta)
        if(max(np.absolute(Beta - Beta_temp))<max_tol):
            break

    res_sum = np.sum(resid**2)
    return res_sum, Beta
