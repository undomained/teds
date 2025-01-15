# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by GM."""
from dataclasses import dataclass
from typing import Self

from pyquaternion import Quaternion
import numpy as np
import numpy.typing as npt


@dataclass
class Navigation:
    """Viewing and solar geometry describing a slice or full orbit."""
    # Timestamps of orbit positions and attitude quaternions, in
    # seconds from beginning of day
    time: npt.NDArray[np.float64]
    # Orbit positions
    orb_pos: npt.NDArray[np.float64]
    # Attitude quaternions
    att_quat: npt.NDArray[Quaternion]
    # Satellite altitude, m
    altitude: npt.NDArray[np.float64]


@dataclass
class Geometry:
    """Viewing and solar geometry describing a slice or full orbit."""
    lat: npt.NDArray[np.float64]
    lon: npt.NDArray[np.float64]
    height: npt.NDArray[np.float64]
    # Solar/viewing azimuth/zenith angles
    saa: npt.NDArray[np.float64]
    sza: npt.NDArray[np.float64]
    vaa: npt.NDArray[np.float64]
    vza: npt.NDArray[np.float64]

    @classmethod
    def from_shape(cls, shape: tuple) -> Self:
        """Initialize all geometry arrays with the same shape."""
        return cls(np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape))
