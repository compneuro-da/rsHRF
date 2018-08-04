import numpy as np
import warnings

warnings.filterwarnings("ignore")

def wgr_get_parameters(hdrf, dt):
    """
    Find Model Parameters
    h - Height
    p - Time to peak (in units of dt where dt = TR/para.T)
    w - Width at half-peak
    """
    param = np.zeros((3, 1))

    if(np.any(hdrf)):
        n = np.fix(np.amax(hdrf.shape) * 0.8)

        p = np.argmax(np.absolute(hdrf[np.arange(0, n, dtype='int')]))
        h = hdrf[p]

        if h > 0:
            v = (hdrf >= (h / 2))
        else:
            v = (hdrf <= (h / 2))
        v = v.astype(int)
        b = np.argmin(np.diff(v))
        v[b + 1:] = 0
        w = np.sum(v)

        cnt = p - 1
        g = hdrf[1:] - hdrf[0:-1]

        while cnt > 0 and np.abs(g[cnt]) < 0.001:
            h = hdrf[cnt - 1]
            p = cnt
            cnt = cnt - 1

        param[0] = h
        param[1] = (p + 1) * dt
        param[2] = w * dt

    else:
        print('.')
    return param.ravel()
