# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Stubs for algorithms."""
import numpy as np
import numpy.typing as npt


def cpp_rt_act(co2_concentration: npt.NDArray[np.float64],
               ch4_concentration: npt.NDArray[np.float64],
               h2o_concentration: npt.NDArray[np.float64],
               co2_xsec: npt.NDArray[np.float64],
               ch4_xsec: npt.NDArray[np.float64],
               h2o_xsec: npt.NDArray[np.float64],
               albedo: npt.NDArray[np.float64],
               mu_sza: npt.NDArray[np.float64],
               mu_vza: npt.NDArray[np.float64],
               sun: npt.NDArray[np.float64],
               rad: npt.NDArray[np.float64]) -> None: ...
