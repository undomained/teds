# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of calibration steps applied to data by the L1A-L1B processor.

Each function should set the new data level accordingly.

"""
from scipy.fft import fft2
from scipy.fft import ifft2
from scipy.interpolate import CubicSpline
from scipy.interpolate import interpn
from tqdm import tqdm
import numpy as np
import numpy.typing as npt

from .binning import bin_data
from .binning import unbin_data
from .types import BinningTable
from .types import CKDNoise
from .types import CKDNonlin
from .types import CKDStray
from .types import CKDSwath
from .types import L1
from .types import ProcLevel


def coadding_and_binning(l1_product: L1,
                         count_table: npt.NDArray[np.int32]) -> None:
    """Undo coadding over time and binning over the detector.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    count_table
        Number of pixels in each bin of a binned frame.

    """
    l1_product.proc_level = ProcLevel.raw
    l1_product.signal /= l1_product.coad_factor
    l1_product.signal /= count_table


def dark_offset(l1_product: L1, offset: npt.NDArray[np.float64]) -> None:
    """Remove offset.

    Parameters
    ----------
      l1_product
        L1 product (signal and detector settings).
      offset
        Detector map of offset [counts].

    """
    l1_product.proc_level = ProcLevel.dark_offset
    l1_product.signal -= offset


def noise(l1_product: L1,
          count_table: npt.NDArray[np.int32],
          ckd: CKDNoise,
          dark_current: npt.NDArray[np.float64],
          artificial_scaling: float) -> None:
    """Determine expected noise in signal.

    The signal is not changed.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Noise CKD consisting of maps of read_noise [counts] and the
        conversion gain [e/counts].

    Returns
    -------
        Detector map of dark current [counts/s].

    """
    l1_product.proc_level = ProcLevel.noise
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = (
        ckd.read_noise**2 + np.abs(l1_product.signal) / ckd.conversion_gain)
    l1_product.noise = artificial_scaling * np.sqrt(
        variance / (l1_product.coad_factor * count_table))


def dark_current(l1_product: L1,
                 dark_current: npt.NDArray[np.float64]) -> None:
    """Remove dark signal.

    If the dark signal does not depend linearly on exposure time, make
    sure the dark current CKD is valid for the used exposure time.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    dark_current
        Detector map of dark current [counts/s].

    """
    l1_product.proc_level = ProcLevel.dark_current
    l1_product.signal -= dark_current * l1_product.exposure_time


def nonlinearity(l1_product: L1,
                 pixel_mask: npt.NDArray[np.bool_],
                 ckd: CKDNonlin) -> None:
    """Remove nonlinearity.

    Only a dependence on signal is taken into account, not on
    intensity, so make sure the non-linearity CKD is valid for the
    used exposure time.

    The model is a linear interpolation of a set of points, so make
    sure the number of points is large enough to produce a smooth
    curve. Input data outside the model range are clipped. The model
    does not depend on pixel.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Nonlinearity CKD consisting of an expected, linear signal
        [counts], and an observed, non-linear signal [counts].

    """
    l1_product.proc_level = ProcLevel.nonlin
    dx = 0.001 * np.min(np.diff(ckd.expected))
    good = ~pixel_mask
    for i_alt in tqdm(range(l1_product.signal.shape[0])):
        l1_product.noise[i_alt, good] *= (
            (np.interp(l1_product.signal[i_alt, good] + dx,
                       ckd.observed,
                       ckd.expected)
             - np.interp(l1_product.signal[i_alt, good] - dx,
                         ckd.observed,
                         ckd.expected))
            / (2 * dx))
        l1_product.signal[i_alt, good] = np.interp(
            l1_product.signal[i_alt, good], ckd.observed, ckd.expected)


def prnu(l1_product: L1,
         pixel_mask: npt.NDArray[np.bool_],
         prnu_qe: npt.NDArray[np.float64]) -> None:
    """Remove PRNU and quantum efficiency.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    prnu_qe
        Detector map of PRNU times quantum efficiency (not the
        correction).

    """
    l1_product.proc_level = ProcLevel.prnu
    good = ~pixel_mask
    for i_alt in tqdm(range(l1_product.signal.shape[0])):
        _prnu = prnu_qe[good]
        l1_product.signal[i_alt, good] /= _prnu
        l1_product.noise[i_alt, good] /= _prnu


def remove_bad_values(n_cols: int,
                      pixel_mask: npt.NDArray[np.bool_],
                      signals: npt.NDArray[np.float64]) -> None:
    """Smooth over bad values in detector images.

    Some algorithms like stray light require all pixels to have a
    normal signal value.

    """
    n_rows = len(pixel_mask) // n_cols
    mask = pixel_mask.reshape((n_rows, n_cols))
    for i_alt in tqdm(range(signals.shape[0])):
        smooth_signal = np.empty((n_rows, n_cols))
        for i_row in range(n_rows):
            good = np.where(~mask[i_row, :])[0]
            signal_row = signals[i_alt, :].reshape((n_rows, n_cols))[i_row, :]
            spline = CubicSpline(good, signal_row[good])
            smooth_signal[i_row, :] = spline(np.arange(n_cols))
        signals[i_alt, :] = smooth_signal.ravel()


def convolve(
        signal: npt.NDArray[np.float64],
        kernel_fft: npt.NDArray[np.complex128]) -> npt.NDArray[np.float64]:
    """Convolve a detector image with kernel.

    Parameters
    ----------
    signal
        Detector image in real space.
    kernel_fft
        Fourier transform of the kernel.

    Returns
    -------
        Convolution result in real space.

    """
    signal = np.pad(signal, ((0, kernel_fft.shape[0]-signal.shape[0]),
                             (0, kernel_fft.shape[1]-signal.shape[1])))
    signal_fft = fft2(signal)
    signal_fft = signal_fft * kernel_fft
    signal_convolved = ifft2(signal_fft).real
    return signal_convolved


def convolve_with_all_kernels(signal: npt.NDArray[np.float64],
                              ckd: CKDStray) -> npt.NDArray[np.float64]:
    """Convolve a detector image with multiple kernels.

    Parameters
    ----------
    signal
        Detector image in real space.
    ckd
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subsignals, and an 'edges' array which
        specifies the location of each subsignal within the original
        signal.

    Returns
    -------
        Convolution result in real space.

    """
    original_shape = signal.shape
    detector_shape = ckd.weights[0, :, :].shape
    signal_convolved = np.zeros(detector_shape)
    for kernel_fft, weights, edges in zip(ckd.kernels_fft,
                                          ckd.weights,
                                          ckd.edges):
        signal_weighted = signal.reshape(detector_shape) * weights
        sub_signal = signal_weighted[edges[0]:edges[1], edges[2]:edges[3]]
        conv_result = convolve(sub_signal, kernel_fft)
        signal_convolved[edges[0]:edges[1], edges[2]:edges[3]] += (
            conv_result[:edges[1]-edges[0], :edges[3]-edges[2]])
    return signal_convolved.reshape(original_shape)


def stray_light(l1_product: L1,
                binning_table: BinningTable,
                ckd: CKDStray,
                van_cittert_steps: int) -> None:
    """Correct a detector image for stray light.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subimages, and an 'edges' array which
        specifies the location of each subimage within the original
        image.
    van_cittert_steps
        Number of deconvolution iterations

    Returns detector images corrected (cleaned) of stray light as part
    of the L1 product.

    """
    l1_product.proc_level = ProcLevel.stray
    if van_cittert_steps == 0:
        return
    n_alt = l1_product.signal.shape[0]
    eta = ckd.eta.ravel()
    for i_alt in tqdm(range(n_alt)):
        signal = unbin_data(binning_table, l1_product.signal[i_alt])
        signal_ideal = signal
        for i_vc in range(van_cittert_steps):
            stray_convolved = convolve_with_all_kernels(signal_ideal, ckd)
            signal_ideal = (signal - stray_convolved) / (1 - eta)
        l1_product.signal[i_alt, :] = bin_data(
            binning_table, signal_ideal, False)


def map_from_detector(l1_product: L1,
                      ckd: CKDSwath,
                      count_table: npt.NDArray[np.int32],
                      wavelengths: npt.NDArray[np.float64],
                      exact_drawing: bool) -> None:
    """Map detector to spectra.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Swath CKD with variables that describe mapping from spectra to
        detector.
    count_table
        Number of pixels in each bin of a binned frame.
    wavelengths
        Wavelength [nm] at a given spectrum and column.

    """
    l1_product.proc_level = ProcLevel.swath
    n_alt = l1_product.signal.shape[0]
    n_cols = ckd.act_map.shape[1]
    n_rows = len(count_table) // n_cols
    n_act, n_wavelengths = ckd.row_map.shape
    l1_product.spectra = np.empty((n_alt, n_act, len(ckd.wavelengths)))
    l1_product.spectra_noise = np.empty(l1_product.spectra.shape)
    row_col_map = np.column_stack((ckd.row_map.ravel(), ckd.col_map.ravel()))
    row_indices = np.zeros(n_rows)
    bin_sizes = count_table.reshape((n_rows, n_cols))[:, 0]
    for i in range(n_rows):
        if i == 0:
            row_indices[i] = 0.5 * bin_sizes[i] - 0.5
        else:
            bin_size_half = 0.5 * (bin_sizes[i-1] - bin_sizes[i])
            row_indices[i] = row_indices[i-1] + bin_sizes[i] + bin_size_half
    col_indices = np.arange(n_cols)
    if not exact_drawing:
        for i_alt in tqdm(range(n_alt)):
            # Compute spectra and noise at intermediate wavelengths
            l1_product.spectra[i_alt, :, :] = interpn(
                (row_indices, col_indices),
                l1_product.signal[i_alt, :].reshape((n_rows, n_cols)),
                row_col_map,
                method='quintic',
                bounds_error=False,
                fill_value=None).reshape((n_act, n_wavelengths))
            l1_product.spectra_noise[i_alt, :, :] = interpn(
                (row_indices, col_indices),
                l1_product.noise[i_alt, :].reshape((n_rows, n_cols)),
                row_col_map,
                method='quintic',
                bounds_error=False,
                fill_value=None).reshape((n_act, n_wavelengths))
        l1_product.wavelengths = ckd.wavelengths
        return
    # When using the exact drawing algorithm, derive spectra directly
    # on the final wavelengths grids (not intermediate).
    l1_product.spectra = np.empty((n_alt, n_act, n_cols))
    l1_product.spectra_noise = np.empty(l1_product.spectra.shape)
    # Regrid row_map and spectra from intermediate wavelengths to
    # wavelengths corresponding to detector columns.
    n_act = ckd.row_map.shape[0]
    row_map = np.empty((n_act, n_cols))
    for i_act in range(n_act):
        row_map[i_act, :] = CubicSpline(
            ckd.wavelengths, ckd.row_map[i_act, :])(wavelengths[i_act, :])
    for i_alt in tqdm(range(n_alt)):
        signal = l1_product.signal[i_alt, :].reshape((n_rows, n_cols))
        noise = l1_product.noise[i_alt, :].reshape((n_rows, n_cols))
        spectra = l1_product.spectra[i_alt, :, :]
        spectra_noise = l1_product.spectra_noise[i_alt, :, :]
        for i_act in range(n_act):
            for i_wave in range(n_cols):
                row_dn = row_map[i_act, i_wave].astype(int)
                row_up = row_dn + 1
                weight = row_map[i_act, i_wave] - row_dn
                spectra[i_act, i_wave] = (
                    weight * signal[row_dn, i_wave]
                    + (1 - weight) * signal[row_up, i_wave])
                spectra_noise[i_act, i_wave] = np.sqrt(
                    (weight * noise[row_dn, i_wave])**2
                    + ((1 - weight) * noise[row_up, i_wave])**2)
    l1_product.wavelengths = wavelengths


def change_wavelength_grid(l1_product: L1,
                           wavelengths_out: npt.NDArray[np.float64]) -> None:
    """Interpolate spectra to new wavelength grid.

    This is required if the instrument model did not draw anything on
    the detector (cal_level was set to swath or higher). Then the
    spectra are still on the intermediate wavelength grid and should
    be interpolated onto the main CKD wavelength grids for the L1B
    product.

    """
    n_alt = l1_product.spectra.shape[0]
    n_act = l1_product.spectra.shape[1]
    new_spectra = np.empty((n_alt, n_act, wavelengths_out.shape[1]))
    new_spectra_noise = np.empty(new_spectra.shape)

    for i_alt in range(n_alt):
        for i_act in range(n_act):
            if len(l1_product.wavelengths.shape) == 2:
                wavelengths = l1_product.wavelengths[i_act, :]
            else:
                wavelengths = l1_product.wavelengths
            spline = CubicSpline(wavelengths,
                                 l1_product.spectra[i_alt, i_act, :])
            new_spectra[i_alt, i_act, :] = spline(wavelengths_out[i_act, :])
            spline = CubicSpline(wavelengths,
                                 l1_product.spectra_noise[i_alt, i_act, :])
            new_spectra_noise[i_alt, i_act, :] = spline(
                wavelengths_out[i_act, :])
    l1_product.spectra = new_spectra
    l1_product.spectra_noise = new_spectra_noise
    l1_product.wavelengths = wavelengths_out


def radiometric(l1_product: L1, rad_corr: npt.NDArray[np.float64]) -> None:
    """Convert counts to spectral photon radiance [nm-1 s-1 sr-1 m-2].

    The quantum efficiency [e/ph] is taken into account by the PRNU
    step.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    rad_corr
        Radiance responsivity correction factor [nm-1 sr-1 m-2].
    exptime
        Exposure time [s] if not already in l1_product.

    """
    l1_product.proc_level = ProcLevel.l1b
    l1_product.spectra *= rad_corr[0, 0] / l1_product.exposure_time
    l1_product.spectra_noise *= rad_corr[0, 0] / l1_product.exposure_time


def bin_l1b(l1_product: L1, bin_spectra: int) -> None:
    """Bin all observables of L1B product.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    bin_spectra
        Binning factor across ACT dimension

    """
    def bin_ACT(arr: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Bin across the ACT dimension.

        Array shape is assumed to be (n_alt, n_act, n_wave), where
        n_wave is the number of wavelengths. Arrays that have a
        different shape need to be reshaped before and after calling
        this function.

        """
        n_alt, n_act, n_wave = arr.shape
        shape_binned = (n_alt, n_act // bin_spectra, n_wave)
        arr_binned = np.zeros(shape_binned)
        for i in range(bin_spectra):
            arr_binned[:, :, :] += arr[:, i::bin_spectra, :]
        arr_binned /= bin_spectra
        return arr_binned

    n_alt, n_act, n_wave = l1_product.spectra.shape
    l1_product.spectra = bin_ACT(l1_product.spectra)
    l1_product.spectra_noise = bin_ACT(l1_product.spectra_noise)
    l1_product.wavelengths = bin_ACT(
        l1_product.wavelengths.reshape(1, n_act, n_wave)).reshape(-1, n_wave)
    l1_product.geometry.lat = bin_ACT(
        l1_product.geometry.lat.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.lon = bin_ACT(
        l1_product.geometry.lon.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.height = bin_ACT(
        l1_product.geometry.height.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.sza = bin_ACT(
        l1_product.geometry.sza.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.saa = bin_ACT(
        l1_product.geometry.saa.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.vza = bin_ACT(
        l1_product.geometry.vza.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
    l1_product.geometry.vaa = bin_ACT(
        l1_product.geometry.vaa.reshape(n_alt, n_act, 1)).reshape(n_alt, -1)
