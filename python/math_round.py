import math
import numpy as np

def np_round_up(x):
    if type(x) is np.ndarray:
        floor_mask = x - np.floor(x) < 0.5
        z = np.where(floor_mask,np.floor(x),np.ceil(x))
        return z.astype(dtype=np.int)
    else:
        if x - math.floor(x) < 0.5:
            return math.floor(x)
        return math.ceil(x)