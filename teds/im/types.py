# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by instrument model."""
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt


@dataclass
class Stokes:
    """Stokes Q and U parameters.

    Radiance (Stokes I) is already part of the L1 class.

    """
    Q: npt.NDArray[np.floating]
    U: npt.NDArray[np.floating]
