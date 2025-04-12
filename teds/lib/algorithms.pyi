# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Stubs for algorithms."""
import numpy as np
import numpy.typing as npt


def cpp_rt_act(co2_concentration: npt.NDArray[np.floating],
               ch4_concentration: npt.NDArray[np.floating],
               h2o_concentration: npt.NDArray[np.floating],
               co2_xsec: npt.NDArray[np.floating],
               ch4_xsec: npt.NDArray[np.floating],
               h2o_xsec: npt.NDArray[np.floating],
               albedo: npt.NDArray[np.floating],
               mu_sza: npt.NDArray[np.floating],
               mu_vza: npt.NDArray[np.floating],
               sun: npt.NDArray[np.floating],
               rad: npt.NDArray[np.floating]) -> None: ...
