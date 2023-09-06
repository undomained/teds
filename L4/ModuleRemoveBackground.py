#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 2023.

@author: Manu Goudar
"""


from scipy.optimize import curve_fit
import numpy as np
from scipy import signal
from .ModuleDataContainer import DataCont


def gauss_linear(x, H0, H1, A, x0, sigma):
    _tmp = -((x - x0) ** 2) / (2 * sigma**2)
    _idx = _tmp > -200
    _tmp1 = np.zeros_like(_tmp)
    _tmp1[_idx] = np.exp(_tmp[_idx])
    return H0 + H1 * x + A * _tmp1


def gauss_fit_linear(x, y, consts=[]):
    _bounds = [(0, -0.001, 0, -5, -10), (np.mean(y), 0.001, np.max(y), 5, 10)]
    try:
        if len(consts) > 0:
            H0 = consts[0]
            H1 = consts[1]
            A = consts[2]
            mean = consts[3]
            sigma = consts[4]
            __p0 = [H0, H1, A, mean, sigma]
            popt, pcov = curve_fit(gauss_linear, x, y, p0=__p0,
                                   check_finite=True, bounds=_bounds,
                                   method="trf", maxfev=15000)
        else:
            H0 = min(y)
            H1 = 0.0001
            A = max(y)
            mean = 0
            sigma = 4
            __p0 = [H0, H1, A, mean, sigma]
            popt, pcov = curve_fit(gauss_linear, x, y, p0=__p0,
                                   check_finite=True, bounds=_bounds,
                                   method="trf", maxfev=15000)
        return popt
    except (RuntimeError, TypeError, NameError):
        return []


def fit_things(x1, y1, info):
    # compute background
    if len(info) == 0:
        H0, H1, A, mean, sigma = gauss_fit_linear(x1, y1)
    else:
        H0, H1, A, mean, sigma = gauss_fit_linear(x1, y1, info)
    y1_fit = gauss_linear(x1, H0, H1, A, mean, sigma)
    return [H0, H1, A, mean, sigma], y1_fit


def pad_arrays(x, y, _pad=0):
    """
    Pad x and y arrays
    """
    # pad y
    _l1 = np.ones(_pad) * y[0]
    _l2 = np.ones(_pad) * y[-1]
    _y = np.concatenate((_l1, y, _l2))
    # pad x
    _l1 = np.arange(-1 * _pad, 0) * 0.5 + x[0]
    _l2 = np.arange(1, _pad + 1) * 0.5 + x[-1]
    _x = np.concatenate((_l1, x, _l2))
    return _x, _y


def get_change_sign(yy, limit=15):
    _tmp = signal.argrelmin(yy, order=5)[0]
    # find the point where things start to increase
    if len(_tmp) == 0:
        return len(yy) - 1
    elif len(_tmp) == 1:
        _n = _tmp[0]
        return _n
    else:
        # check if it is within the limit and if it is switch it to next
        # if the value increased by really small amount (1%) or
        # if the values go down in next 4 indices w.r.t value at index _n
        if _tmp[0] < limit:
            return _tmp[1]
        else:
            return _tmp[0]


def get_center_and_cutoff_index(ydata, limit=10):
    """
    Center the data and get index of changing sign
    """
    nos = np.int_(len(ydata) / 2)
    _tmp = np.where(ydata == np.max(ydata[nos - limit: nos + limit]))[0]
    nos = _tmp[len(_tmp) // 2]
    yy1 = np.flip(ydata[: nos + 1])
    yy2 = ydata[nos:]
    _idx1 = get_change_sign(yy1)
    _idx2 = get_change_sign(yy2)
    return nos, [nos - _idx1, nos + _idx2 + 1]


def emission_indices(x, y):
    """
    Get indices for data to be used for computing emissions
    """
    _ix = np.argwhere(x == 0)[0][0]
    fst = np.where(y[:_ix] < 0)[0]
    lst = np.where(y[_ix:] < 0)[0]
    if len(fst) > 0:
        i1 = fst[-1]
    else:
        i1 = 0
    if len(lst) > 0:
        i2 = lst[0] + _ix
    else:
        i2 = len(y)
    _idx = np.zeros_like(y, dtype=np.bool_)
    _idx[i1+1: i2] = True
    return _idx


def remove_background(initdata, diff=0.15):
    """Compute enhancement by removing background along a transect.

    Parameters
    ----------
    initdata : classmethod
        Initial data of a transaction line.
    diff : type
        Description of parameter `diff` (the default is 0.15).

    Returns
    -------
    None
    """
    # center and cut off the plume
    c_id, idx = get_center_and_cutoff_index(initdata.co)

    # cut off data
    data = DataCont()
    data.__setattr__("cutoff_co", initdata.co[idx[0]:idx[1]])
    data.__setattr__("cutoff_line", initdata.line[idx[0]:idx[1]] -
                     initdata.line[c_id])
    data.__setattr__("cutoff_coords", initdata.coords[idx[0]:idx[1]])

    # Get cut off indices and pad the array for background computation
    x11 = data.cutoff_line.copy()/1000   # convert the mts to km
    y11 = data.cutoff_co.copy()
    x1, y1 = pad_arrays(x11, y11, 20)

    # If the difference between two sides is not high then continue
    if abs(y1[0] - y1[-1]) / np.max(y1) < diff:
        _params = np.zeros((6))
        _params[0] = (y1[0] - y1[-1]) / np.max(y1)
        # fit a gaussian curve and compute background removed CO
        _params[1:], yfit1 = fit_things(x1, y1, [])
        data.__setattr__("gaussfit_x", x1)
        data.__setattr__("gaussfit_co", yfit1)
        ynew_back_removed = y11 - (_params[1] + _params[2] * x11)
        data.__setattr__("gaussfit_params", _params)
        data.__setattr__("back_removed", ynew_back_removed)

        eid = emission_indices(x11, data.back_removed)
        if len(ynew_back_removed[eid]) > 0:
            data.__setattr__('filtered_index', eid)
            data.__setattr__('backgroundremovalflag', True)
        else:
            data.__setattr__('backgroundremovalflag', False)
    else:
        data.__setattr__("backgroundremovalflag", False)
    return data
