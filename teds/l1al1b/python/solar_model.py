# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Wrapper for the solar model.

"""
from pyquaternion import Quaternion
import numpy as np
import numpy.typing as npt

try:
    from .geolocation import solar_model as c_solar_model
except ModuleNotFoundError:
    def c_solar_model(tai_seconds: np.uint,
                      tai_second_fraction: np.floating,
                      q_ecef_j2000: npt.NDArray[np.floating]) -> None:
        pass


def solar_model(tai_seconds: np.uint,
                tai_second_fraction: np.floating) -> Quaternion:
    q_ecef_j2000 = np.empty(4)
    c_solar_model(tai_seconds, tai_second_fraction, q_ecef_j2000)
    return Quaternion(q_ecef_j2000[3], *q_ecef_j2000[:3]).normalised
