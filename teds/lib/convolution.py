# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions related to numerical convolution."""
from scipy.signal import convolve
import math
import numba
import numpy as np
import numpy.typing as npt


class KernelGauss:
    """Kernel for generalized Gaussian convolution.

    The kernel is stored as a matrix and convolution is performed as
    matrix-vector multiplication.

    """
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
        # Number of rows in each kernel
        n_rows = wave_out.shape[0]
        # Matrix representing the convolution kernels
        self.kernel = np.empty((len(wave_out), len(wave_in)), dtype=np.float64)
        # Ranges of non-zero values for each kernel row
        self.i_limits = np.empty((n_rows, 2), dtype=np.int32)
        self._gen_kernel(fwhm,
                         shape,
                         cutoff,
                         wave_in,
                         wave_out,
                         self.i_limits,
                         self.kernel)

    @staticmethod
    @numba.njit
    def _gen_kernel(fwhm: float,
                    shape: float,
                    cutoff: float,
                    wave_in: npt.NDArray[np.float64],
                    wave_out: npt.NDArray[np.float64],
                    i_limits: npt.NDArray[np.int32],
                    kernel: npt.NDArray[np.float64]) -> None:
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

    @staticmethod
    @numba.njit
    def _convolve(i_limits: npt.NDArray[np.int32],
                  kernel: npt.NDArray[np.float64],
                  array: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Convolve the kernel with an array."""
        conv = np.empty(kernel.shape[0])
        for i_row in range(kernel.shape[0]):
            _lim = i_limits[i_row, :]
            conv[i_row] = np.dot(kernel[i_row, _lim[0]:_lim[1]],
                                 array[_lim[0]:_lim[1]])
            # conv[i_row] = np.dot(
            #     kernel[i_row, i_limits[i_row, 0]:i_limits[i_row, 1]],
            #     array[i_limits[i_row, 0]:i_limits[i_row, 1]])
        return conv

    def convolve(self,
                 array: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Convolve the kernel with an array."""
        return KernelGauss._convolve(self.i_limits, self.kernel, array)


class KernelGauss2D:
    """Kernel for generalized 2D Gaussian convolution.

    The kernel is stored as a small matrix and convolution is
    performed with FFTs.

    """
    def __init__(self,
                 fwhm_x: float,
                 fwhm_y: float,
                 shape_x: float,
                 shape_y: float) -> None:
        dim_scaling = 6
        dim_x = max(int(dim_scaling * fwhm_x), 1)
        dim_y = max(int(dim_scaling * fwhm_y), 1)
        # Make sure dimension is an odd number so the kernel can
        # always be centered.
        dim_x += 1 - dim_x % 2
        dim_y += 1 - dim_y % 2
        x0 = dim_x // 2
        y0 = dim_y // 2
        x = 2**(-(2 * np.abs(np.arange(dim_x) - x0) / fwhm_x)**(shape_x))
        y = 2**(-(2 * np.abs(np.arange(dim_y) - y0) / fwhm_y)**(shape_y))
        self.kernel = np.outer(y, x)
        self.kernel /= self.kernel.sum()

    def convolve(self,
                 data: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return convolve(data, self.kernel, mode='same')
