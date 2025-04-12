# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by GM."""
from dataclasses import dataclass
from pyquaternion import Quaternion
from typing import Self
import numpy as np
import numpy.typing as npt


@dataclass
class Navigation:
    """Viewing and solar geometry describing a slice or full orbit."""
    # Timestamps of orbit positions and attitude quaternions, in
    # seconds from beginning of day
    time: npt.NDArray[np.floating]
    # Orbit positions
    orb_pos: npt.NDArray[np.floating]
    # Attitude quaternions
    att_quat: npt.NDArray[Quaternion]
    # Satellite altitude, m
    altitude: npt.NDArray[np.floating]

    @classmethod
    def from_shape(cls, shape: tuple) -> Self:
        """Initialize all geometry arrays with the same shape."""
        return cls(np.zeros(shape),
                   np.zeros(shape),
                   np.zeros(shape, dtype=Quaternion),
                   np.zeros(shape))


@dataclass
class Geometry:
    """Viewing and solar geometry describing a slice or full orbit."""
    lat: npt.NDArray[np.floating]
    lon: npt.NDArray[np.floating]
    height: npt.NDArray[np.floating]
    # Solar/viewing azimuth/zenith angles
    sza: npt.NDArray[np.floating]
    saa: npt.NDArray[np.floating]
    vza: npt.NDArray[np.floating]
    vaa: npt.NDArray[np.floating]

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

    def deg2rad(self) -> None:
        self.lat = np.deg2rad(self.lat)
        self.lon = np.deg2rad(self.lon)
        self.sza = np.deg2rad(self.sza)
        self.saa = np.deg2rad(self.saa)
        self.vza = np.deg2rad(self.vza)
        self.vaa = np.deg2rad(self.vaa)

    def rad2deg(self) -> None:
        self.lat = np.rad2deg(self.lat)
        self.lon = np.rad2deg(self.lon)
        self.sza = np.rad2deg(self.sza)
        self.saa = np.rad2deg(self.saa)
        self.vza = np.rad2deg(self.vza)
        self.vaa = np.rad2deg(self.vaa)
