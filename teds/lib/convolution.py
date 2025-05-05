# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions related to numerical convolution."""
from netCDF4 import Dataset
from scipy.interpolate import CubicSpline
from scipy.signal import convolve
from typing import Self
import math
import numba
import numpy as np
import numpy.typing as npt


class Kernel:
    """Class for convolutions with a kernel.

    The init and other functions below refer to an ISRF and wavelength
    grids because that is a typical use case for convolutions. But
    this class can also be used with other types of information (in
    general, some input/output grids and kernel data). The
    ISRF/wavelength terminology is easier to follow.

    Input and output grids are assumed to be uniform.

    """
    def __init__(self,
                 wavelength_diffs: npt.NDArray[np.floating],
                 isrf: npt.NDArray[np.floating],
                 wavelengths_in: npt.NDArray[np.floating],
                 wavelengths_out: npt.NDArray[np.floating],
                 wave_cutoff: float = 0.58) -> None:
        """Precompute kernel values that will be used for convolution.

        Parameters
        ----------
        wavelength_diffs
            Wavelengths at which the ISRF (kernel) data is given
            (centered around 0)
        isrf
            ISRF data corresponding to wavelength_diffs. Can be
            different from input/output grids used for convolution.
        wavelengths_in
            Wavelength grid of unconvolved data
        wavelengths_out
            Wavelength grid of convolved data
        wave_cutoff
            Extent of the ISRF, in nm. If the kernel is (well)
            localized then a small value of this parameter can be used
            to reduced computational cost. The dimension needs to be
            the same as wavelength_diffs so not necessarily nm.

        """
        self.wavelength_diffs = wavelength_diffs
        self.isrf_data = isrf
        # Convolution wavelength range
        self.wavelengths_in = wavelengths_in
        self.wavelengths_out = wavelengths_out
        self.wave_step = self.wavelengths_in[1] - self.wavelengths_in[0]
        # This many values of the kernel need to be precomputed and
        # stored. Outside the range the kernel is 0.
        self.n_vals_half = int(wave_cutoff // self.wave_step)
        self.n_vals = 2 * self.n_vals_half - 1
        self.regenerate()

    def regenerate(self, x0: float = 0) -> None:
        # Interpolate input data onto the wavelength range that is
        # relevant for convolution.
        isrf_spline = CubicSpline(self.wavelength_diffs - x0, self.isrf_data)
        self.isrf_values = isrf_spline(-self.n_vals_half * self.wave_step
                                       + range(self.n_vals) * self.wave_step)
        norm = self.isrf_values.sum()
        self.isrf_values = self.isrf_values / norm
        self.isrf_values_der = isrf_spline.derivative()(
            -self.n_vals_half * self.wave_step
            + range(self.n_vals) * self.wave_step)
        self.isrf_values_der /= norm

    @classmethod
    def from_file(cls,
                  filename: str,
                  wavelengths_in: npt.NDArray[np.floating],
                  wavelengths_out: npt.NDArray[np.floating]) -> Self:
        """Read ISRF from file."""
        nc = Dataset(filename)
        return cls(nc['wavelength'][:].data,
                   nc['isrf'][:].data,
                   wavelengths_in,
                   wavelengths_out)

    @classmethod
    def from_gauss(cls,
                   wavelengths_in: npt.NDArray[np.floating],
                   wavelengths_out: npt.NDArray[np.floating],
                   fwhm: float,
                   shape: float = 2.0,
                   wave_cutoff: float = 0.7) -> Self:
        """Generate kernel from generalized Gaussian parameters.

        The generalized Gaussian analytical form is used to construct
        the wavelength range and the ISRF values which are then passed
        to the normal constructor. It's basically a convenience
        function to generate the ISRF data. Otherwise, this class
        works the same as with tabulated data.

        Parameters
        ----------
        wavelengths_in
            Wavelength grid of unconvolved data
        wavelengths_out
            Wavelength grid of convolved data
        fwhm
            FWHM of the generalized Gaussian
        shape
            Shape parameter of the generalized Gaussian (2 means
            normal Gaussian). Lower values lead to stronger wings and
            larger values to a more blocky shape.
        wave_cutoff
            Kernel extent

        Returns
        -------
            Instance of the kernel class

        """
        step = wavelengths_in[1] - wavelengths_in[0]
        wavelength_diffs = np.arange(-wave_cutoff, wave_cutoff, step)
        isrf = 2**(-(2 * np.abs(wavelength_diffs) / fwhm)**(shape))
        return cls(wavelength_diffs, isrf, wavelengths_in, wavelengths_out)

    @staticmethod
    @numba.njit
    def _convolve(wavelengths_in: npt.NDArray[np.floating],
                  wavelengths_out: npt.NDArray[np.floating],
                  wave_step: float,
                  n_vals_half: npt.NDArray[np.floating],
                  n_vals: int,
                  isrf_values: npt.NDArray[np.floating],
                  data: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Convolve array.

        This needs to be static method to use Numba decorator. The
        arguments are not actually variable by the calling function.

        Parameters
        ----------
        wavelengths_in
            Wavelengths of unconvolved data
        wavelengths_out
            Wavelengths of convolved data
        wave_step
            Step size of wavelengths_in
        n_vals_half
            Half the number of ISRF evaluations required for convolution
        n_vals
            Full number of ISRF evaluations required for convolution
        isrf_values
            ISRF values tabulated for n_vals elements in steps of wave_step
        data
            Data to be convolved

        Returns
        -------
            Convolved array

        """
        conv = np.empty(len(wavelengths_out))
        # Find index of first element in the unconvolved wavelength
        # array corresponding to a target (convolved data)
        # wavelength. Do this for all wavelengths at once.
        first_in_wavelength = wavelengths_out - n_vals_half * wave_step
        first_in_idx = np.searchsorted(wavelengths_in, first_in_wavelength) - 1
        for i_wave in range(len(wavelengths_out)):
            idx_beg = first_in_idx[i_wave]
            conv[i_wave] = np.dot(isrf_values, data[idx_beg:idx_beg+n_vals])
        return conv

    def convolve(self,
                 data: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Convolve a data array with the kernel.

        Parameter
        ---------
        data
            Data to be convolved. The associated grid is
            wavelengths_in as specified in the constructor and cannot
            be changed anymore.

        Returns
        -------
            Convolved array. The assocated grid is wavelengths_out.

        """
        return Kernel._convolve(self.wavelengths_in,
                                self.wavelengths_out,
                                self.wave_step,
                                self.n_vals_half,
                                self.n_vals,
                                self.isrf_values,
                                data)

    def convolve_der(
            self,
            data: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Convolve with the derivative of the kernel."""
        return Kernel._convolve(self.wavelengths_in,
                                self.wavelengths_out,
                                self.wave_step,
                                self.n_vals_half,
                                self.n_vals,
                                self.isrf_values_der,
                                data)


class KernelGauss:
    """Kernel for generalized Gaussian convolution.

    The kernel is stored as a matrix and convolution is performed as
    matrix-vector multiplication.

    """
    def __init__(self,
                 wave_in: npt.NDArray[np.floating],
                 wave_out: npt.NDArray[np.floating],
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
        self.kernel = np.empty((len(wave_out), len(wave_in)),
                               dtype=np.floating)
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
                    wave_in: npt.NDArray[np.floating],
                    wave_out: npt.NDArray[np.floating],
                    i_limits: npt.NDArray[np.int32],
                    kernel: npt.NDArray[np.floating]) -> None:
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
                  kernel: npt.NDArray[np.floating],
                  array: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Convolve the kernel with an array."""
        conv = np.empty(kernel.shape[0])
        for i_row in range(kernel.shape[0]):
            _lim = i_limits[i_row, :]
            conv[i_row] = np.dot(kernel[i_row, _lim[0]:_lim[1]],
                                 array[_lim[0]:_lim[1]])
        return conv

    def convolve(self,
                 array: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
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
                 data: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        return convolve(data, self.kernel, mode='same')
