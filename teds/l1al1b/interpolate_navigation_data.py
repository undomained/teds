# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Interpolate platform navigation data to detector timestamps."""
from pyquaternion import Quaternion
from scipy.interpolate import CubicSpline
import numpy as np
import numpy.typing as npt

from .types import Navigation


def interpolate_navigation_data(
        navigation: Navigation,
        timestamps: npt.NDArray[np.floating]) -> Navigation:
    """Interpolate navigation data from orbit times to detector image
    times.

    The returned navigation data will be part of the L1 product class.

    Parameters
    ----------
    navigation
        Platform navigation data.
    timestamps
        Detector timestamps onto which to interpolate navigation data,
        in seconds from beginning of day.

    Returns
    -------
        Navigation data interpolated to detector timestamps.

    """
    n_alt = len(timestamps)
    orb_pos = np.empty((n_alt, 3))
    att_quat = np.empty(n_alt, dtype=Quaternion)
    altitude = np.zeros((n_alt,))
    # Interpolate orbit positions
    for i_dir in range(3):
        s = CubicSpline(navigation.time, navigation.orb_pos[:, i_dir])
        orb_pos[:, i_dir] = s(timestamps)
    # Interpolate quaternions
    indices = np.searchsorted(navigation.time.astype(np.float64), timestamps)
    for i_alt in range(n_alt):
        idx_lo = indices[i_alt] - 1
        idx_hi = indices[i_alt]
        idx_delta = ((timestamps[i_alt] - navigation.time[idx_lo])
                     / (navigation.time[idx_hi] - navigation.time[idx_lo]))
        q0 = navigation.att_quat[idx_lo]
        q1 = navigation.att_quat[idx_hi]
        att_quat[i_alt] = Quaternion.slerp(q0, q1, amount=idx_delta)
    return Navigation(timestamps, orb_pos, att_quat, altitude)
