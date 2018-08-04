import numpy as np
import warnings

warnings.filterwarnings("ignore")

def knee_pt(y, x=None):
    x_was_none = False
    use_absolute_dev_p = True
    res_x = np.nan
    idx_of_result = np.nan

    if type(y) is not np.ndarray:
        print('knee_pt: y must be a numpy 1D vector')
        return res_x, idx_of_result
    else:
        if y.ndim >= 2:
            print('knee_pt: y must be 1 dimensional')
            return res_x, idx_of_result
        if np.size(y) == 0:
            print('knee_pt: y can not be an empty vector')
            return res_x, idx_of_result
        else:
            if x is None:
                x_was_none = True
                x = np.arange(1, np.amax(y.shape) + 1, dtype=np.int)
            if x.shape != y.shape:
                print('knee_pt: y and x must have the same dimensions')
                return res_x, idx_of_result
            if y.size < 3:
                print('knee_pt: y must be at least 3 elements long')
                return res_x, idx_of_result
            if np.all(np.diff(x) >= 0) and (not x_was_none):
                idx = np.argsort(x)
                y = np.sort(y)
                x = np.sort(x)
            else:
                idx = np.arange(0, np.amax(x.shape))
            sigma_xy = np.cumsum(np.multiply(x, y), axis=0)
            sigma_x = np.cumsum(x, axis=0)
            sigma_y = np.cumsum(y, axis=0)
            sigma_xx = np.cumsum(np.multiply(x, x), axis=0)
            n = np.arange(1, np.amax(y.shape) + 1).conj().T
            det = np.multiply(n, sigma_xx) - np.multiply(sigma_x, sigma_x)
            mfwd = (np.multiply(n, sigma_xy) -
                    np.multiply(sigma_x, sigma_y)) / det
            bfwd = -1 * ((np.multiply(sigma_x, sigma_xy) -
                          np.multiply(sigma_xx, sigma_y)) / det)

            sigma_xy = np.cumsum(np.multiply(x[::-1], y[::-1]), axis=0)
            sigma_x = np.cumsum(x[::-1], axis=0)
            sigma_y = np.cumsum(y[::-1], axis=0)
            sigma_xx = np.cumsum(np.multiply(x[::-1], x[::-1]), axis=0)
            n = np.arange(1, np.amax(y.shape) + 1).conj().T
            det = np.multiply(n, sigma_xx) - np.multiply(sigma_x, sigma_x)
            mbck = ((np.multiply(n, sigma_xy) -
                     np.multiply(sigma_x, sigma_y)) / det)[::-1]
            bbck = (-1 *
                    ((np.multiply(sigma_x, sigma_xy) -
                      np.multiply(sigma_xx, sigma_y)) / det))[::-1]

            error_curve = np.full(y.shape, np.nan)
            for breakpt in range(1, np.amax((y - 1).shape)):
                delsfwd = (np.multiply(mfwd[breakpt], x[0:breakpt + 1]) +
                           bfwd[breakpt]) - y[0:breakpt + 1]
                delsbck = (np.multiply(mbck[breakpt], x[breakpt:]) +
                           bbck[breakpt]) - y[breakpt:]
                if use_absolute_dev_p:
                    error_curve[breakpt] = \
                        np.sum(np.abs(delsfwd)) + np.sum(np.abs(delsbck))
                else:
                    error_curve[breakpt] = \
                        np.sqrt(np.sum(np.multiply(delsfwd, delsfwd))) + \
                        np.sqrt(np.sum(np.multiply(delsbck, delsbck)))
            try:
                loc = np.nanargmin(error_curve)
            except ValueError as e:
                loc = 0
            res_x = x[loc]
            idx_of_result = idx[loc]
    return res_x, idx_of_result
