# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Stubs for geolocation."""
import numpy as np
import numpy.typing as npt


def solar_model(tai_seconds: np.uint,
                tai_second_fraction: np.float64,
                q_ecef_j2000: npt.NDArray[np.float64]) -> None: ...


def geolocate(dem_filename: str,
              los: npt.NDArray[np.float64],
              tai_seconds: npt.NDArray[np.uint],
              tai_subsec: npt.NDArray[np.float64],
              orb_pos: npt.NDArray[np.float64],
              att_quat: npt.NDArray[np.float64],
              lat: npt.NDArray[np.float64],
              lon: npt.NDArray[np.float64],
              height: npt.NDArray[np.float64],
              vza: npt.NDArray[np.float64],
              vaa: npt.NDArray[np.float64],
              sza: npt.NDArray[np.float64],
              saa: npt.NDArray[np.float64]) -> None: ...
