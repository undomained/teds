# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of calibration steps applied to data by the L1A-L1B processor.

Each function should set the new data level accordingly.

"""
from scipy.fft import fft2
from scipy.fft import ifft2
import numpy as np
import numpy.typing as npt


from .binning import unbin_data
from .types import CKDNoise
from .types import CKDNonlin
from .types import CKDStray
from .types import L1
from .types import ProcLevel


def remove_coadding_and_binning(l1_product: L1,
                                bin_indices: npt.NDArray[np.int32],
                                count_table: npt.NDArray[np.int32],
                                method: str) -> None:
    """Undo coadding over time and binning over the detector.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    bin_indices
        Bin index of each pixel in an unbinned frame.
    count_table
        Number of pixels in each bin of a binned frame.
    method
        Method for unbinning the data, but 'none' keeps the data
        binned, divided by the binning factors and bins the CKD. Other
        options are 'nearest' to split each bin value equally between
        all corresponding pixels (keeping the sum over all data the
        same), 'linear' or 'cubic' to put each ratio of bin value and
        binning factor at the centroid of the corresponding pixels and
        interpolate.

    """
    # C++ code does this while reading data
    l1_product['signal']
    l1_product['signal'] /= l1_product['coad_factors'][..., None]
    # Signal corresponds now with exposure time, but
    # l1_product['coad_factors'] is not adjusted to remember which
    # coaddition factors were used.
    if method.lower() == 'none':
        # Binning unchanged but sum changed to mean for later
        # non-linearity correction.
        l1_product['signal'] /= count_table
        # Binning factors, not noise, C++ code writes 0s
        l1_product['noise'] = np.tile(count_table,
                                      (len(l1_product['signal']), 1))
    else:
        l1_product['signal'] = unbin_data(
            l1_product['signal'], bin_indices, count_table, method=method)
        # Binning table id is kept in sync with the data
        l1_product['binning_table_ids'] = np.zeros_like(
            l1_product['binning_table_ids'])
        # Original binning factors, C++ code writes 0s
        l1_product['noise'] = unbin_data(
            np.tile(count_table**2, (len(l1_product['signal']), 1)),
            bin_indices,
            count_table)
    l1_product['proc_level'] = ProcLevel.raw


def remove_offset(l1_product: L1, offset: npt.NDArray[np.float64]) -> None:
    """Remove offset.

    Parameters
    ----------
      l1_product
        L1 product (signal and detector settings).
      offset
        Detector map of offset [counts].

    """
    # Assuming bad pixels are NaN already
    l1_product['signal'] -= offset.ravel()
    # Uncertainty of offset not implemented and would mess with the
    # binning factors in RAW data. No change of processing level.


def determine_noise(l1_product: L1,
                    ckd: CKDNoise,
                    pixel_mask: npt.NDArray[bool],
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
    orig_bin_factors = l1_product['noise']  # quirk of RAW data
    # The signal is split in contributions from light and dark in case
    # the dark signal is negative.
    dark_signal = dark_current.ravel() * l1_product['exptimes'][..., None]
    photoresponse_signal = l1_product['signal'] - dark_signal
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = ckd['read_noise'].ravel()**2 + (
        photoresponse_signal + dark_signal)/ckd['conversion_gain'].ravel()
    l1_product['noise'] = np.sqrt(
        np.clip(variance, 0, None)/l1_product['coad_factors'][..., None])
    if l1_product['binning_table_ids'][0] > 0:
        l1_product['noise'] /= np.sqrt(orig_bin_factors)
    # C++ code flags pixels with a nonpositive noise, but there are
    # two problems with that:
    # * in general, bad pixels should be flagged in the CKD
    # * in this case, a nonpositive noise is more likely a model problem
    l1_product['noise'][variance <= 0] = np.nan
    for i_alt in range(variance.shape[0]):
        pixel_mask[variance[i_alt, :].reshape(pixel_mask.shape) <= 0] = True


def remove_dark_signal(l1_product: L1,
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
    # Assuming bad pixels are NaN already
    l1_product['signal'] -= (
        dark_current.ravel() * l1_product['exptimes'][..., None])
    # uncertainty of dark signal not implemented
    l1_product['proc_level'] = ProcLevel.dark


def remove_nonlinearity(l1_product: L1, ckd: CKDNonlin) -> None:
    """Remove non-linearity.

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
    # C++ code gives noise that is orders of magnitude too large
    dx = 0.001 * np.min(np.diff(ckd['expected']))
    l1_product['noise'] *= (
        (np.interp(l1_product['signal'] + dx,
                   ckd['observed'],
                   ckd['expected'])
         - np.interp(l1_product['signal'] - dx,
                     ckd['observed'],
                     ckd['expected']))
        / (2 * dx))
    # C++ code gives extremely large values instead of clipping
    l1_product['signal'] = np.interp(
        l1_product['signal'], ckd['observed'], ckd['expected'])
    l1_product['proc_level'] = ProcLevel.nonlin


def remove_prnu(l1_product: L1, prnu_qe: npt.NDArray[np.float64]) -> None:
    """Remove PRNU and quantum efficiency.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    prnu_qe
        Detector map of PRNU times quantum efficiency (not the
        correction).

    """
    # Assuming bad pixels are NaN already
    l1_product['signal'] /= prnu_qe.ravel()
    l1_product['noise'] /= prnu_qe.ravel()
    l1_product['proc_level'] = ProcLevel.prnu


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


def remove_bad_values(pixel_mask: npt.NDArray[bool],
                      image: npt.NDArray[np.float64]) -> None:
    mask = pixel_mask.ravel()
    for i_alt in range(image.shape[0]):
        for i in range(image.shape[1]):
            if mask[i]:
                n_good = 0
                image[i_alt, i] = 0.0
                if i > 0 and not mask[i - 1]:
                    image[i_alt, i] += image[i_alt, i-1]
                    n_good += 1
                if i < image.shape[1] - 1 and not mask[i + 1]:
                    image[i_alt, i] += image[i_alt, i+1]
                    n_good += 1
                image[i_alt, i] /= max(1, n_good)


def stray_light(l1_product: L1, ckd: CKDStray, van_cittert_steps: int) -> None:
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

    Returns
    -------
        Detector image corrected (cleaned) of stray light.

    """
    if van_cittert_steps == 0:
        return
    eta = ckd['eta'].reshape(l1_product['signal'][0, :].shape[0])
    for i_alt in range(l1_product['signal'].shape[0]):
        image = l1_product['signal'][i_alt]
        image_ideal = l1_product['signal'][i_alt]
        for i_vc in range(van_cittert_steps):
            stray_convolved = convolve_with_all_kernels(image_ideal, ckd)
            image_ideal = (image - stray_convolved) / (1 - eta)
        l1_product['signal'][i_alt, :] = image_ideal


def map_from_detector(l1_products: L1,
                      wavelengths: npt.NDArray[np.float64],
                      spectrum_rows: npt.NDArray[np.float64],
                      spectrum_width: float,
                      pixel_mask: npt.NDArray[bool],
                      bin_indices: npt.NDArray[np.int32]) -> None:
    """Map detector to spectra.

    Spectra are collected from curves on the detector with a fixed
    width along columns. This is not the exact inverse of
    `map_to_detector`.

    Parameters
    ----------
    l1_products
        L1 products (signal and detector settings).
    wavelengths
        Wavelength [nm] at a given spectrum and column.
    spectrum_rows
        Central row index (float) of a given spectrum and at each
        column, 0 is halfway first row.
    spectrum_width
        Number of pixels along a column centred around spectrum_rows
        belonging to one spectrum, at least 1.0.
    pixel_mask
        Detector map of pixel flags, False good, True
        bad.
    bin_indices
        Bin index of each pixel in an unbinned frame.

    """
    # spectrum_indices: indices of associated detector pixels for a
    # given spectrum and column, referring to the detector as a 1D
    # array
    # spectrum_weights: relative weights of the associated detector
    # pixels, normalized or NaN
    n_spectra, n_cols = spectrum_rows.shape
    n_indices = int(np.ceil(spectrum_width))+1
    # Pixel indices referring to the detector as a 1D array
    spectrum_indices = np.empty((n_spectra, n_cols, n_indices), dtype=int)
    spectrum_weights = np.ones((n_spectra, n_cols, n_indices))
    for s in range(n_spectra):
        for c in range(n_cols):
            # 0 is outer edge first row
            row_first = spectrum_rows[s, c] + 0.5 - spectrum_width/2
            row_last = spectrum_rows[s, c] + 0.5 + spectrum_width/2
            spectrum_indices[s, c] = (int(row_first)
                                      + np.arange(n_indices))*n_cols+c
            spectrum_weights[s, c, -1] = 0
            spectrum_weights[s, c, 0] = np.ceil(row_first)-row_first
            spectrum_weights[s, c, int(row_last)-int(row_first)] = (
                row_last-np.floor(row_last))
    if l1_products['binning_table_ids'][0] > 0:
        spectrum_indices = bin_indices.ravel()[spectrum_indices]
    # C++ code takes also bad pixels into account that are flagged
    # during processing, but currently that only happens during the
    # noise determination (see there).
    spectrum_weights *= (
        1.0 - np.take(pixel_mask, spectrum_indices).astype(float))
    weight_sums = np.sum(spectrum_weights, axis=2)
    weight_sums[weight_sums < 0.01] = np.nan  # C++ code uses == 0
    spectrum_weights /= weight_sums[..., None]
    l1_products['spectra'] = np.array(
        [np.sum(spectrum_weights*np.take(signal, spectrum_indices), axis=2)
         for signal in l1_products['signal']])
    l1_products['spectra_noise'] = np.sqrt(
        np.array([np.sum((
            spectrum_weights*np.take(noise, spectrum_indices))**2,
                         axis=2) for noise in l1_products['noise']]))
    l1_products['wavelengths'] = wavelengths
    l1_products['proc_level'] = ProcLevel.swath


def convert_to_radiance(l1_products: L1,
                        rad_corr: npt.NDArray[np.float64],
                        exptime: npt.NDArray[np.float64]) -> None:
    """Convert counts to spectral photon radiance [nm-1 s-1 sr-1 m-2].

    The quantum efficiency [e/ph] is taken into account by the PRNU
    step.

    Parameters
    ----------
    l1_products
        L1 products (signal and detector settings).
    rad_corr
        Radiance responsivity correction factor [nm-1 sr-1 m-2].
    exptime
        Exposure time [s] if not already in l1_products.

    """
    # should never be the case
    if np.isnan(l1_products['exptimes']).all():
        l1_products['exptimes'] = np.full_like(l1_products['exptimes'],
                                               exptime)
    if np.min(l1_products['exptimes']) <= 0:
        raise SystemExit("ERROR: exposure time cannot be "
                         f"{np.min(l1_products['exptimes'])} s")
    # mean_gain = np.mean(conversion_gain)
    mean_gain = 1  # in C++ code effectively 1
    l1_products['spectra'] *= (
        rad_corr*mean_gain/l1_products['exptimes'][:, None, None])
    l1_products['spectra_noise'] *= (
        rad_corr*mean_gain/l1_products['exptimes'][:, None, None])
    l1_products['proc_level'] = ProcLevel.l1b
