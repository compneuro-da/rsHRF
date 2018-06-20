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
