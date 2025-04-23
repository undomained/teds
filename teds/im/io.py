# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions for reading and writing instrument model data."""
from netCDF4 import Dataset
import numpy as np
import numpy.typing as npt

from .types import Stokes


def read_pol(sgm_filename: str,
             mueller_filename: str,
             alt_beg: int,
             alt_end: int | None,) -> tuple[Stokes | None,
                                            npt.NDArray[np.floating]]:
    """Read Stokes Q and U parameters and Mueller matrix elements."""
    if alt_end is not None:
        alt_end += 1
    nc = Dataset(sgm_filename)
    if 'q' in nc.variables:
        stokes = Stokes(nc['q'][alt_beg:alt_end].data,
                        nc['u'][alt_beg:alt_end].data)
    else:
        stokes = None
    if mueller_filename:
        nc = Dataset(mueller_filename)
        mueller = nc['mueller'][:, alt_beg:alt_end].data
    else:
        mueller = np.empty((0))
    return stokes, mueller
