# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Stubs for geolocation."""
import numpy as np
import numpy.typing as npt


def solar_model(tai_seconds: np.uint,
                tai_second_fraction: np.floating,
                q_ecef_j2000: npt.NDArray[np.floating]) -> None: ...


def geolocate(dem_filename: str,
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
              saa: npt.NDArray[np.floating]) -> None: ...
