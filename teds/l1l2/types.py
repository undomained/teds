# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by L2 processor."""
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt


@dataclass
class RefProfiles:
    def __init__(self) -> None:
        self.gases: dict[str, npt.NDArray[np.floating]] = {}
        self.initial: dict[str, float] = {}


class L2:
    """Level 2 data product with diagnostics properties."""
    def __init__(self,
                 n_alt: int,
                 n_act: int,
                 n_wave: int,
                 n_lay: int,
                 n_albedo: int,
                 gas_names: list[str]) -> None:
        # Retrieval diagnostics
        self.chi2 = np.empty((n_alt, n_act))
        self.converged = np.full((n_alt, n_act), False, dtype=np.bool_)
        self.iterations = np.empty((n_alt, n_act), dtype=np.int32)

        # Albedo coefficients, dry air mixing ratios (e.g. XCO2),
        # proxy mixing ratios, and their precisions, accuracies, and
        # other related quantities. Each of those variables is a
        # dictionary with Numpy arrays. This way gases can easily be
        # fetched by their name.
        self.albedo0 = np.empty((n_alt, n_act))
        self.mixing_ratios: dict[str, npt.NDArray[np.floating]] = {}
        self.precisions = {}
        self.gains = {}
        self.col_avg_kernels = {}
        self.proxys: dict[str, npt.NDArray[np.floating]] = {}
        self.proxy_precisions = {}
        for gas in gas_names:
            self.mixing_ratios[gas] = np.empty((n_alt, n_act))
            self.precisions[gas] = np.empty((n_alt, n_act))
            self.gains[gas] = np.empty((n_alt, n_act, n_wave))
            self.col_avg_kernels[gas] = np.empty((n_alt, n_act, n_lay))
            if gas in ('CO2', 'CH4'):
                self.proxys[gas] = np.empty((n_alt, n_act))
                self.proxy_precisions[gas] = np.empty((n_alt, n_act))

        # Unused for now
        self.spec_shift = np.zeros((n_alt, n_act))
        self.spec_squeeze = np.zeros((n_alt, n_act))
