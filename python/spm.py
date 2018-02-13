import nibabel as nib
import numpy as np
from scipy.special import gammaln

def spm_vol(input_nii_file):
    """
    Get header information for images
    """
    v = nib.load(input_nii_file)
    return v

def spm_read_vols(mapped_image_volume):
    """
    Read in entire image volumes
    """
    data = mapped_image_volume.get_data()
    data = data.flatten(order='F')
    return data

def spm_orth(X,OPT='pad'):
    """
    Recursive Gram-Schmidt orthogonalisation of basis functions
    @X - matrix
    @OPT - 'norm' - for Euclidean normalisation
           'pad'  - for zero padding of null space (default)
    """
    def gs_cofficient(v1, v2):
        return np.dot(v2, v1) / np.dot(v1, v1)

    def multiply(cofficient, v):
        return map((lambda x : x * cofficient), v)

    def proj(v1, v2):
        return multiply(gs_cofficient(v1, v2) , v1)

    def gs(X, row_vecs=True, norm = True):
        if not row_vecs:
            X = X.T
        Y = X[0:1,:].copy()
        for i in range(1, X.shape[0]):
            proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
            Y = np.vstack((Y, X[i,:] - proj.sum(0)))
        if norm:
            Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
        if row_vecs:
            return Y
        else:
            return Y.T

    if OPT == 'norm':
        return gs(X,row_vecs=False,norm=True)
    elif OPT == 'pad':
        return gs(X,row_vecs=False,norm=False)
    else:
        return X

def spm_hrf(RT, P=None, fMRI_T=16):
    """
    @RT - scan repeat time
    @P  - parameters of the response function (two gamma functions)

    defaults  (seconds)
    %	P[0] - Delay of Response (relative to onset)	    6
    %	P[1] - Delay of Undershoot (relative to onset)     16
    %	P[2] - Dispersion of Response			            1
    %	P[3] - Dispersion of Undershoot			            1
    %	P[4] - Ratio of Response to Undershoot		        6
    %	P[5] - Onset (seconds)				                0
    %	P[6] - Length of Kernel (seconds)		           32

    hrf  - hemodynamic response function
    P    - parameters of the response function
    """
    p = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    if P is not None:
        p[0:len(P)] = P
    _spm_Gpdf = lambda x, h, l: np.exp(h * np.log(l) + (h - 1) * np.log(x) - (l * x) - gammaln(h))
    # modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / float(fMRI_T)
    u = np.arange(0, int(p[6] / dt + 1)) - p[5] / dt
    with np.errstate(divide='ignore'):  # Known division-by-zero
        hrf = _spm_Gpdf(u, p[0] / p[2], dt / p[2]) - _spm_Gpdf(u, p[1] / p[3],
                                                               dt / p[3]) / p[4]
    idx = np.arange(0, int((p[6] / RT) + 1)) * fMRI_T
    hrf = hrf[idx]
    hrf = np.nan_to_num(hrf)
    hrf = hrf / np.sum(hrf)
    return hrf
