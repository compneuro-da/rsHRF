import numpy as np
from spm import *
from scipy.sparse import lil_matrix
from scipy import stats

def wgr_hrf_estimation_canon(dat,xBF,length,N,bf,temporal_mask):
    """
    Estimate HRF
    """
    thr = xBF['thr']
    u0 = wgr_BOLD_event_vector(N,dat,thr,temporal_mask)
    u = np.append(u0.toarray(),np.zeros((xBF['T']-1,N)),axis=0)
    u = np.reshape(u,(1,-1),order = 'F')
    beta, lag = wgr_hrf_fit(dat,length,xBF,u,N,bf)
    beta_hrf = beta
    beta_hrf = np.append(beta_hrf,lag)
    return beta_hrf, u0

def wgr_BOLD_event_vector(N,matrix,thr,temporal_mask):
    """
    Detect BOLD event.
    event > thr & event < 3.1
    """
    data = lil_matrix((1,N))
    matrix = matrix[:,np.newaxis]
    if 0 in np.array(temporal_mask).shape:
        matrix = stats.zscore(matrix)
        for t in range(3,N-1):
            if ( matrix[t-1,0] > thr and matrix[t-1,0] < 3.1 and np.all(matrix[t-3:t-1,0] < matrix[t-1,0]) and np.all(matrix[t-1,0] > matrix[t:t+2,0]) ):
                data[0,t-1] = 1
    else:
        datm = np.mean(matrix[temporal_mask])
        datstd = np.std(matrix[temporal_mask])
        datstd[datstd==0] = 1
        matrix = np.divide((matrix - datm),datstd)
        for t in range(3,N-1):
            if temporal_mask[t]:
                if ( matrix[t-1,0] > thr and matrix[t-1,0] < 3.1 and np.all(matrix[t-3:t-1,0] < matrix[t-1,0]) and np.all(matrix[t-1,0] > matrix[t:t+2,0]) ):
                    data[0,t-1] = 1
    return data

def wgr_hrf_fit(dat,length,xBF,u,N,bf):
    """
    @u    - BOLD event vector (microtime).
    @nlag - time lag from neural event to BOLD event
    """
    lag = xBF['lag']
    AR_lag = xBF['AR_lag']
    nlag = len(lag)
    erm = np.zeros((1,nlag))
    beta = np.zeros((bf.shape[1]+1,nlag))
    for i in range(nlag):
        u_lag = np.append(u[0,lag[i]:],np.zeros((1,lag[i]))).T
        erm[0,i], beta[:,i] = wgr_glm_estimation(dat,u_lag,bf,xBF['T'],xBF['T0'],AR_lag)

    idx = np.argmin(erm)
    return beta[:,idx], lag[idx]

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
