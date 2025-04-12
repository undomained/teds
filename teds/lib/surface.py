# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tools to simulate surface reflection"""
import numpy as np
import numpy.typing as npt


class Surface:
    """Collection of methods to calculate surface properties"""
    def __init__(self, wave: npt.NDArray[np.float64]) -> None:
        wave_max = np.amax(wave)
        wave_min = np.amin(wave)
        wave_mean = np.mean(wave)
        spec = (wave - wave_mean) / (wave_max - wave_min)
        self.spec = spec
        self.alb = np.zeros_like(wave)

    def get_albedo_poly(self, alb_coeff: list[float]) -> None:
        """Get albedo.

        Parameters
        ----------
        spec
            Spectral array (can be any spectral quantity)
        alb
            List of albedo coefficients

        Returns
        -------
            Albedo polynomial dependent on wavelength.

        """
        albedo = np.zeros(len(self.spec))
        for i in range(len(alb_coeff)):
            albedo += alb_coeff[i] * self.spec**i
        self.alb = albedo
