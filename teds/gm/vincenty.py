# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Vincenty algorithm for computing the distance between two geodetic
points.

"""
import numpy as np
import numpy.typing as npt


def vincenty(lat1: npt.NDArray[np.float64],
             lat2: npt.NDArray[np.float64],
             lon1: npt.NDArray[np.float64],
             lon2: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Vincenty algorithm.

    Parameters
    ----------
    lat1
        Latitude of source point (radians)
    lat2
        Latitude of destination point (radians)
    lon1
        Longitude of source point (radians)
    lon2
        Longitude of destination point (radians)

    Returns
    -------
        Distance, in meters, between the two points along the Earth's
        curvature.

    """
    a = 6378137.0
    b = 6356752.314245179499
    f = 3.3528106647474804385e-3
    U1 = np.arctan((1 - f) * np.tan(lat1))
    U2 = np.arctan((1 - f) * np.tan(lat2))
    sU1 = np.sin(U1)
    sU2 = np.sin(U2)
    cU1 = np.cos(U1)
    cU2 = np.cos(U2)
    L = lon2 - lon1
    cos_alpha2 = 0
    sigma = 0
    sin_sigma = 0
    cos_sigma = 0
    cos_2sig = 0
    _lambda = L
    _lambda_prev = _lambda
    max_iter = 10
    for i_iter in range(max_iter):
        slam = np.sin(_lambda)
        clam = np.cos(_lambda)
        cos_U2_t1 = cU2 * slam
        cos_U2_t2 = cU1 * sU2 - sU1 * cU2 * clam
        sin_sigma = np.sqrt(cos_U2_t1 * cos_U2_t1 + cos_U2_t2 * cos_U2_t2)
        cos_sigma = sU1 * sU2 + cU1 * cU2 * clam
        sigma = np.arctan2(sin_sigma, cos_sigma)
        sin_alpha = np.where(np.abs(sin_sigma < 1e-30),
                             cU1,
                             cU1 * cU2 * slam / sin_sigma)
        cos_alpha2 = 1 - sin_alpha * sin_alpha
        cos_2sig = cos_sigma - 2 * sU1 * sU2 / cos_alpha2
        C = f / 16 * cos_alpha2 * (4.0 + f * (4 - 3 * cos_alpha2))
        _lambda = (L + (1 - C) * f * sin_alpha
                   * (sigma + C * sin_sigma
                      * (cos_2sig + C * cos_sigma * (-1 + 2 * cos_2sig))))
        if (abs(_lambda - _lambda_prev).max() < 1e-12):
            break
        _lambda_prev = _lambda
    u2 = cos_alpha2 * (a * a - b * b) / (b * b)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    cos_2sig2 = cos_2sig * cos_2sig
    D_sigma = (B * sin_sigma
               * (cos_2sig
                  + 1.0 / 4 * B
                  * (cos_sigma * (-1 + 2 * cos_2sig2)
                     - B / 6 * cos_2sig
                     * (-3 + 4 * sin_sigma * sin_sigma)
                     * (-3 + 4 * cos_2sig2))))
    # This is a mypy bug. The type is correct.
    return b * A * (sigma - D_sigma)  # type: ignore
