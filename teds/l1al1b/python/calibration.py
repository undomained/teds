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
    l1_product['proc_level'] = ProcLevel.raw
    l1_product['image'] /= l1_product['coad_factors'][..., None]
    l1_product['image'] /= count_table


def dark_offset(l1_product: L1, offset: npt.NDArray[np.float64]) -> None:
    """Remove offset.

    Parameters
    ----------
      l1_product
        L1 product (signal and detector settings).
      offset
        Detector map of offset [counts].

    """
    l1_product['proc_level'] = ProcLevel.dark_offset
    l1_product['image'] -= offset


def noise(l1_product: L1,
          count_table: npt.NDArray[np.int32],
          ckd: CKDNoise,
          dark_current: npt.NDArray[np.float64]) -> None:
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
    l1_product['proc_level'] = ProcLevel.noise
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = (ckd['read_noise']**2
                + np.abs(l1_product['image']) / ckd['conversion_gain'])
    l1_product['noise'] = np.sqrt(
        variance / (l1_product['coad_factors'][..., None] * count_table))


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
    l1_product['proc_level'] = ProcLevel.dark_current
    l1_product['image'] -= dark_current * l1_product['exptimes'][..., None]


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
    l1_product['proc_level'] = ProcLevel.nonlin
    dx = 0.001 * np.min(np.diff(ckd['expected']))
    good = ~pixel_mask
    for i_alt in tqdm(range(l1_product['image'].shape[0])):
        l1_product['noise'][i_alt, good] *= (
            (np.interp(l1_product['image'][i_alt, good] + dx,
                       ckd['observed'],
                       ckd['expected'])
             - np.interp(l1_product['image'][i_alt, good] - dx,
                         ckd['observed'],
                         ckd['expected']))
            / (2 * dx))
        l1_product['image'][i_alt, good] = np.interp(
            l1_product['image'][i_alt, good],
            ckd['observed'],
            ckd['expected'])


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
    l1_product['proc_level'] = ProcLevel.prnu
    good = ~pixel_mask
    for i_alt in tqdm(range(l1_product['image'].shape[0])):
        _prnu = prnu_qe[good]
        l1_product['image'][i_alt, good] /= _prnu
        l1_product['noise'][i_alt, good] /= _prnu


def convolve(
        image: npt.NDArray[np.float64],
        kernel_fft: npt.NDArray[np.complex128]) -> npt.NDArray[np.float64]:
    """Convolve an image with kernel.

    Parameters
    ----------
    image
        Image in real space.
    kernel_fft
        Fourier transform of the kernel.

    Returns
    -------
        Convolution result in real space.

    """
    image = np.pad(image, ((0, kernel_fft.shape[0]-image.shape[0]),
                           (0, kernel_fft.shape[1]-image.shape[1])))
    image_fft = fft2(image)
    image_fft = image_fft * kernel_fft
    image_convolved = ifft2(image_fft).real
    return image_convolved


def convolve_with_all_kernels(image: npt.NDArray[np.float64],
                              ckd: CKDStray) -> npt.NDArray[np.float64]:
    """Convolve an image with multiple kernels.

    Parameters
    ----------
    image
        Image in real space.
    ckd
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subimages, and an 'edges' array which
        specifies the location of each subimage within the original
        image.

    Returns
    -------
        Convolution result in real space.

    """
    original_shape = image.shape
    detector_shape = ckd['weights'][0, :, :].shape
    image_convolved = np.zeros(detector_shape)
    for kernel_fft, weights, edges in zip(ckd['kernels_fft'],
                                          ckd['weights'],
                                          ckd['edges']):
        image_weighted = image.reshape(detector_shape) * weights
        sub_image = image_weighted[edges[0]:edges[1], edges[2]:edges[3]]
        conv_result = convolve(sub_image, kernel_fft)
        image_convolved[edges[0]:edges[1], edges[2]:edges[3]] += (
            conv_result[:edges[1]-edges[0], :edges[3]-edges[2]])
    return image_convolved.reshape(original_shape)


def remove_bad_values(n_cols: int,
                      pixel_mask: npt.NDArray[np.bool_],
                      images: npt.NDArray[np.float64]) -> None:
    """Smooth over bad values in detector images.

    Some algorithms like stray light require all pixels to have a
    normal signal value.

    """
    n_rows = len(pixel_mask) // n_cols
    mask = pixel_mask.reshape((n_rows, n_cols))
    for i_alt in tqdm(range(images.shape[0])):
        smooth_image = np.empty((n_rows, n_cols))
        for i_row in range(n_rows):
            good = np.where(~mask[i_row, :])[0]
            image_row = images[i_alt, :].reshape((n_rows, n_cols))[i_row, :]
            spline = CubicSpline(good, image_row[good])
            smooth_image[i_row, :] = spline(np.arange(n_cols))
        images[i_alt, :] = smooth_image.ravel()


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
    l1_product['proc_level'] = ProcLevel.stray
    if van_cittert_steps == 0:
        return
    n_alt = l1_product['image'].shape[0]
    eta = ckd['eta'].ravel()
    for i_alt in tqdm(range(n_alt)):
        image = unbin_data(binning_table, l1_product['image'][i_alt])
        image_ideal = image
        for i_vc in range(van_cittert_steps):
            stray_convolved = convolve_with_all_kernels(image_ideal, ckd)
            image_ideal = (image - stray_convolved) / (1 - eta)
        l1_product['image'][i_alt, :] = bin_data(binning_table, image_ideal)


def map_from_detector(l1_product: L1,
                      ckd: CKDSwath,
                      count_table: npt.NDArray[np.int32],
                      wavelengths: npt.NDArray[np.float64]) -> None:
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
    l1_product['proc_level'] = ProcLevel.swath
    n_alt = l1_product['image'].shape[0]
    n_cols = ckd['act_map'].shape[1]
    n_rows = len(count_table) // n_cols
    n_act, n_wavelengths = ckd['row_map'].shape
    l1_product['spectra'] = np.empty((n_alt, n_act, len(ckd['wavelengths'])))
    l1_product['spectra_noise'] = np.empty(l1_product['spectra'].shape)
    row_col_map = np.column_stack((ckd['row_map'].ravel(),
                                   ckd['col_map'].ravel()))
    row_indices = np.zeros(n_rows)
    bin_sizes = count_table.reshape((n_rows, n_cols))[:, 0]
    for i in range(n_rows):
        if i == 0:
            row_indices[i] = 0.5 * bin_sizes[i] - 0.5
        else:
            bin_size_half = 0.5 * (bin_sizes[i-1] - bin_sizes[i])
            row_indices[i] = row_indices[i-1] + bin_sizes[i] + bin_size_half
    col_indices = np.arange(n_cols)
    for i_alt in tqdm(range(n_alt)):
        # Compute spectra and noise at intermediate wavelengths
        l1_product['spectra'][i_alt, :, :] = interpn(
            (row_indices, col_indices),
            l1_product['image'][i_alt, :].reshape((n_rows, n_cols)),
            row_col_map,
            method='quintic',
            bounds_error=False,
            fill_value=None).reshape((n_act, n_wavelengths))
        l1_product['spectra_noise'][i_alt, :, :] = interpn(
            (row_indices, col_indices),
            l1_product['noise'][i_alt, :].reshape((n_rows, n_cols)),
            row_col_map,
            method='quintic',
            bounds_error=False,
            fill_value=None).reshape((n_act, n_wavelengths))
    l1_product['wavelengths'] = ckd['wavelengths']


def change_wavelength_grid(l1_product: L1,
                           wavelengths_out: npt.NDArray[np.float64]) -> None:
    """Interpolate spectra to new wavelength grid.

    This is required if the instrument model did not draw anything on
    the detector (cal_level was set to swath or higher). Then the
    spectra are still on the intermediate wavelength grid and should
    be interpolated onto the main CKD wavelength grids for the L1B
    product.

    """
    n_alt = l1_product['spectra'].shape[0]
    n_act = l1_product['spectra'].shape[1]
    new_spectra = np.empty((n_alt, n_act, wavelengths_out.shape[1]))
    new_spectra_noise = np.empty(new_spectra.shape)
    for i_alt in range(n_alt):
        for i_act in range(n_act):
            spline = CubicSpline(l1_product['wavelengths'],
                                 l1_product['spectra'][i_alt, i_act, :])
            new_spectra[i_alt, i_act, :] = spline(wavelengths_out[i_act, :])
            spline = CubicSpline(l1_product['wavelengths'],
                                 l1_product['spectra_noise'][i_alt, i_act, :])
            new_spectra_noise[i_alt, i_act, :] = spline(
                wavelengths_out[i_act, :])
        l1_product['spectra'] = new_spectra
        l1_product['spectra_noise'] = new_spectra_noise
        l1_product['wavelengths'] = wavelengths_out


def radiometric(l1_product: L1,
                rad_corr: npt.NDArray[np.float64],
                exptime: npt.NDArray[np.float64]) -> None:
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
    l1_product['proc_level'] = ProcLevel.l1b
    l1_product['spectra'] *= (
        rad_corr[0, 0] / l1_product['exptimes'][:, None, None])
    l1_product['spectra_noise'] *= (
        rad_corr[0, 0] / l1_product['exptimes'][:, None, None])
