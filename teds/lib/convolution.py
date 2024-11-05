# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions related to numerical convolution."""
import math
import numpy as np
import numpy.typing as npt
import numba


class Kernel:
    """Kernel for generalized Gaussian convolution"""
    def __init__(self,
                 wave_in: npt.NDArray[np.float64],
                 wave_out: npt.NDArray[np.float64],
                 fwhm: float,
                 shape: float = 2.0,
                 cutoff: float = 1.5) -> None:
        """Construct the kernel and optimize summation limits.

        Parameters
        ----------
        wave_in
            Source wavelength grid (must be evenly spaced)
        wave_out
            Target wavelength (can have arbitrary spacing)
        fwhm
            Full width at half maximum of the Gaussian
        shape
            Shape parameter of the generalized normal distribution.
            Value 2 for Gauss (default), towards 1 for stronger wings
            and large values for a more blocky shape.
        cutoff
            Source and target wavelength differences, in multiples of
            FWHM, that are considered for summing. The aim is to
            reduce the computational cost by considering the
            localization of the Gaussian.

        """
        # Number of different kernels
        n_kernels = wave_out.shape[0]
        # Number of rows in each kernel
        n_rows = wave_out.shape[1]
        # Matrix representing the convolution kernels
        self.kernels = np.empty((n_kernels, n_rows, len(wave_in)),
                                dtype=np.float64)
        # Ranges of non-zero values for each kernel row
        self.i_limits = np.empty((n_kernels, n_rows, 2), dtype=np.int32)
        for i_kernel in range(n_kernels):
            self._gen_kernel(fwhm,
                             shape,
                             cutoff,
                             wave_in,
                             wave_out[i_kernel],
                             self.i_limits[i_kernel],
                             self.kernels[i_kernel])
        n_nonzeros = self.i_limits[:, :, 1] - self.i_limits[:, :, 0]
        max_nonzero = n_nonzeros.flatten().max()
        self.kernels = self.kernels[:, :, :max_nonzero]

    @staticmethod
    @numba.njit
    def _gen_kernel(fwhm, shape, cutoff, wave_in, wave_out, i_limits, kernel):
        do_gaussian = abs(shape - 2.0) < 1e-30
        sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        sigma_2inv = 1 / (math.sqrt(2.0) * sigma)
        for i_row in range(kernel.shape[0]):
            diff = wave_in - wave_out[i_row]
            i_limits[i_row, :] = (
                np.searchsorted(diff, -cutoff * fwhm),
                np.searchsorted(diff, cutoff * fwhm))
            _lim = i_limits[i_row, :]
            if do_gaussian:
                kernel[i_row, _lim[0]:_lim[1]] = np.exp(
                    -(diff[_lim[0]:_lim[1]] * sigma_2inv)**shape)
            else:
                kernel[i_row, _lim[0]:_lim[1]] = 2**(
                    -(2 * np.abs(diff[_lim[0]:_lim[1]]) / fwhm)**(shape))
            # Normalize to 1
            kernel[i_row, _lim[0]:_lim[1]] /= np.sum(
                kernel[i_row, _lim[0]:_lim[1]])
            n_nonzero = _lim[1]-_lim[0]
            kernel[i_row, :n_nonzero] = kernel[i_row, _lim[0]:_lim[1]]

    def convolve(self,
                 array: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Convolve the kernel with an array."""
        n_kernels, n_rows, n_cols = self.kernels.shape
        conv = np.empty((n_kernels, n_rows))
        for i_kernel in range(n_kernels):
            for i_row in range(n_rows):
                _lim = self.i_limits[i_kernel, i_row, :]
                conv[i_kernel, i_row] = np.dot(
                    self.kernels[i_kernel, i_row, :_lim[1]-_lim[0]],
                    array[i_kernel, _lim[0]:_lim[1]])
        return conv
