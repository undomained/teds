#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 2023.

@author: Manu Goudar


Usage: One of these
    precision = level4precision_lsq(conc, noise, cdfweight)
    precision = level4precision_lsq_threshold(conc, noise, threshold)
"""

import numpy as np
from .ModulePlume import cdfthresholdbasedplume


def lsq_precision(error_variance, ydata):
    """Computes the level 4 precision by LSQ

    Args:
        error_variance (vector float): Variance of the noise.
        ydata (vector float): Actual data (not retrieved).

    Returns:
        precision (float): Level 4 precision from LSQ
    """
    sy = np.diag(error_variance)
    syinv = np.linalg.inv(sy)
    K = np.reshape(ydata, (ydata.size, 1))
    return np.sqrt(np.linalg.inv(np.dot(K.T, np.dot(syinv, K))))[0][0]


def lsq_emission(error_variance, ydata, noisy):
    """Computes the level 4 precision by LSQ

    Args:
        error_variance (vector float): Variance of the noise.
        ydata (vector float): Actual data (not retrieved).
        NOISY (vector float): Level 2 retrieved data.

    Returns:
        emission (float): emission from LSQ
    """
    sy = np.diag(error_variance)
    syinv = np.linalg.inv(sy)
    K = np.reshape(ydata, (ydata.size, 1))
    sx = np.linalg.inv(np.dot(K.T, np.dot(syinv, K)))
    gain = np.dot(sx, np.dot(K.T, syinv))
    return (np.dot(gain, noisy)).ravel()[0]


def level4precision_lsq(conc, level2_precision, cdfweight):
    """Computes level 4 precision based on % of top pixels.

    Args:
        conc (matrix float): True concentrations from simulated data
                             (not satellite retreival).
        noise (matrix float): Sigma noise per pixel
        cdfweight (float): Enhanced pixels above this value. For example
                           value of 0.8 corresponds to values above 80%.

    Returns:
        precision (float): Level 4 precision from LSQ
    """
    # compute plume
    plume = cdfthresholdbasedplume(conc, cdfweight)
    # compute variance
    variance = level2_precision[plume]**2
    # compute precision
    return lsq_precision(variance, conc[plume])


def level4precision_lsq_threshold(conc, noise, threshold):
    """Computes level 4 precision based on threshold.

    Args:
        conc (matrix float): True concentrations.
        noise (matrix float): Sigma noise per pixel
        threshold (float): Same units as conc.

    Returns:
        precision (float): Level 4 precision from LSQ
    """
    # compute plume
    plume = conc > threshold
    # compute the variance
    variance = noise[plume]**2
    # return precicion
    return lsq_precision(variance, conc[plume])


def emissionprecision(actualconc, level2, level2_precision, cdfweight):


    """Computes level 4 precision based on % of top pixels.

    Args:
        actualconc (matrix float): True concentrations from simulated data
                                   (not satellite retreival).
        level2 (matrix float): True concentrations from simulated data
                               (not satellite retreival).
        lvel2precision (matrix float): Sigma noise per pixel
        cdfweight (float): Enhanced pixels above this value. For example
                           value of 0.8 corresponds to values above 80%.

    Returns:
        emission (float): Emission estimate
        precision (float): Level 4 precision from LSQ
    """
    # compute plume
    plume = cdfthresholdbasedplume(level2, cdfweight)
    # compute variance
    variance = level2_precision[plume]**2
    # compute precision
    precision = lsq_precision(variance, actualconc[plume])
    emission = lsq_emission(variance, actualconc[plume], level2[plume])
    return emission, precision

