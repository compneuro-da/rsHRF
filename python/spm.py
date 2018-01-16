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
