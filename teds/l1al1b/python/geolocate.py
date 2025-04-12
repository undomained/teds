# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Wrapper for L1B processor geolocation routines."""
import numpy as np
import numpy.typing as npt

from .types import Geometry
from .types import L1

try:
    from .geolocation import geolocate as c_geolocate
except ModuleNotFoundError:
    def c_geolocate(dem_filename: str,
                    los: npt.NDArray[np.floating],
                    tai_seconds: npt.NDArray[np.uint],
                    tai_subsec: npt.NDArray[np.floating],
                    orb_pos: npt.NDArray[np.floating],
                    att_quat: npt.NDArray[np.floating],
                    lat: npt.NDArray[np.floating],
                    lon: npt.NDArray[np.floating],
                    height: npt.NDArray[np.floating],
                    vza: npt.NDArray[np.floating],
                    vaa: npt.NDArray[np.floating],
                    sza: npt.NDArray[np.floating],
                    saa: npt.NDArray[np.floating]) -> None:
        pass


def geolocate(
        l1: L1, los: npt.NDArray[np.floating], dem_filename: str) -> Geometry:
    """C

    Parameters
    ----------
    l1
        L1 product. It does not contain any science data (detector
        images or spectra) but has detector image timestamps and
        navigation data required for geolocation.
    los
        Line-of-sight vectors
    dem_filename
        Path to a NetCDF containing the digital elevation model

    Returns
    -------

    """
    n_alt = l1.navigation.orb_pos.shape[0]
    n_act = los.shape[0]
    geometry = Geometry.from_shape((n_alt, n_act))
    # Need to convert Pyquaternion objects into Numpy arrays first
    att_quat = np.empty((n_alt, 4))
    for i_alt in range(n_alt):
        att_quat[i_alt, :] = np.roll(
            l1.navigation.att_quat[i_alt].elements, -1)
    c_geolocate(dem_filename,
                los,
                l1.tai_seconds,
                l1.tai_subsec,
                l1.navigation.orb_pos,
                att_quat,
                geometry.lat,
                geometry.lon,
                geometry.height,
                geometry.vza,
                geometry.vaa,
                geometry.sza,
                geometry.saa)
    return geometry
