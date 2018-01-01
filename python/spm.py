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
