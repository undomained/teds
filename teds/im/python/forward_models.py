# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of processes applied to data by the instrument model.

Each process undoes or "uncalibrates" a step that is performed by the
L1A-L1B processor, gradually bringing the data level from L1B to L1A
or anywhere in between.

"""
from scipy.interpolate import interpn
from tqdm import tqdm
import numpy as np
import numpy.typing as npt

from teds.l1al1b.python.calibration import convolve_with_all_kernels
from teds.l1al1b.python.types import BinningTable
from teds.l1al1b.python.types import CKDNoise
from teds.l1al1b.python.types import CKDNonlin
from teds.l1al1b.python.types import CKDStray
from teds.l1al1b.python.types import CKDSwath
from teds.l1al1b.python.types import ProcLevel
from teds.l1al1b.python.types import L1
from teds.lib import convolution


def apply_isrf(l1_product: L1,
               wavelengths_out: npt.NDArray[np.float64],
               convolve: bool,
               fwhm: float,
               shape: float,
               in_memory: bool) -> None:
    """Convolve spectra with the ISRF.

    The ISRF has a fixed shape as a function of wavelength and does not depend
    on pixel.

    For debugging or experimentation, the ISRF convolution can be
    turned off in which case this function interpolates onto a target
    grid without convolving.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    wavelengths_out
        New wavelength grid [nm].
    convolve
        Convolve with ISRF (True) or only interpolate (False)
    fwhm
        Full width at half-maximum [nm] of ISRF.
    shape
        ISRF shape parameter. Value 2 for Gauss (default), towards 1
        for stronger wings and large values for a more blocky shape.
    in_memory
        Whether to store the entire ISRF in memory or compute it per
        spectrum

    """
    l1_product['proc_level'] = ProcLevel.l1b
    wavelengths_in = l1_product['wavelengths']  # SGM data, so 1D
    n_alt = l1_product['spectra'].shape[0]
    n_act = l1_product['spectra'].shape[1]
    conv = np.empty((n_alt, n_act, len(wavelengths_out)))
    wavelengths_out_tiled = np.tile(wavelengths_out, n_act).reshape(n_act, -1)
    if convolve:
        # Convolve spectra with ISRF assuming wavelengths_in is monotonically
        # increasing or decreasing (checked while
        # reading). Extrapolated values are close to zero assuming no
        # bad values in the spectra. In C++ code fixed to Gauss.
        if in_memory:
            kernel = convolution.Kernel(
                wavelengths_in, wavelengths_out_tiled, fwhm, shape, in_memory)
            for i_alt in tqdm(range(n_alt)):
                conv[i_alt, :, :] = kernel.convolve(
                    l1_product['spectra'][i_alt, :, :])
        else:
            for i_alt in tqdm(range(n_alt), unit=' ALT'):
                for i_act in range(n_act):
                    kernel = convolution.Kernel(
                        wavelengths_in, wavelengths_out, fwhm, shape)
                    conv[i_alt, i_act:i_act+1, :] = kernel.convolve(
                        l1_product['spectra'][i_alt, i_act:i_act+1, :])
        l1_product['spectra'] = conv
    else:
        for i_alt in tqdm(range(n_alt), unit=' ALT'):
            for i_act in range(n_act):
                conv[i_alt, i_act, :] = np.interp(
                    wavelengths_out,
                    wavelengths_in,
                    l1_product['spectra'][i_alt, i_act, :])
    l1_product['spectra'] = conv
    l1_product['wavelengths'] = wavelengths_out_tiled


def radiometric(l1_product: L1, rad_corr: npt.NDArray[np.float64]) -> None:
    """Convert from spectral photon radiance [nm-1 s-1 sr-1 m-2] to
    counts.

    The quantum efficiency [e/ph] is taken into account by the PRNU
    step.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    rad_corr
        Radiance responsivity correction factor [nm-1 sr-1 m-2].

    """
    l1_product['proc_level'] = ProcLevel.swath
    for i_alt in tqdm(range(l1_product['spectra'].shape[0])):
        l1_product['spectra'][i_alt, :, :] *= (
            l1_product['exptimes'][i_alt] / (rad_corr[0, 0]))


def map_to_detector(l1_product: L1, ckd: CKDSwath) -> None:
    """Map spectra to detector.

    Spectra are mapped to infinitely thin curves on the detector (instead of
    simply to rows due to spatial distortion or keystone) and then interpolated
    along columns (not along constant wavelengths) at integer row coordinates,
    using cubic interpolation and linear extrapolation.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    spectrum_rows
        Central row index (float) of a given spectrum and at each
        column.
    n_detrows
        Number of detector rows.

    """
    l1_product['proc_level'] = ProcLevel.stray
    n_alt = l1_product['spectra'].shape[0]
    n_rows, n_cols = ckd['act_map'].shape
    l1_product['image'] = np.empty((n_alt, n_rows * n_cols))
    act_wavelength_map = np.column_stack((ckd['act_map'].ravel(),
                                          ckd['wavelength_map'].ravel()))
    for i_alt in tqdm(range(n_alt)):
        l1_product['image'][i_alt, :] = interpn(
            (ckd['act_angles'], ckd['wavelengths']),
            l1_product['spectra'][i_alt, :, :],
            act_wavelength_map,
            method='quintic',
            bounds_error=False,
            fill_value=None).reshape(n_rows * n_cols)


def stray_light(l1_product: L1, ckd: CKDStray) -> None:
    """Add stray light to the image.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subimages, and an 'edges' array which
        specifies the location of each subimage within the original
        image.

    """
    l1_product['proc_level'] = ProcLevel.prnu
    n_alt = l1_product['image'].shape[0]
    eta = ckd['eta'].ravel()
    for i_alt in tqdm(range(n_alt)):
        image = l1_product['image'][i_alt]
        stray_convolved = convolve_with_all_kernels(image, ckd)
        image_convolved = (1 - eta) * image + stray_convolved
        l1_product['image'][i_alt, :] = image_convolved


def prnu(l1_product: L1, prnu_qe: npt.NDArray[np.float64]) -> None:
    """Incorporate PRNU and quantum efficiency.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    prnu_qe
        Detector map of PRNU times quantum efficiency (not the
        correction).

    """
    l1_product['proc_level'] = ProcLevel.nonlin
    l1_product['image'] *= prnu_qe


def nonlinearity(l1_product: L1, ckd: CKDNonlin) -> None:
    """Incorporate non-linearity.

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
    l1_product['proc_level'] = ProcLevel.dark_current
    l1_product['image'] = np.interp(
        l1_product['image'], ckd['expected'], ckd['observed'])


def dark_current(l1_product: L1,
                 dark_current: npt.NDArray[np.float64]) -> None:
    """Incorporate dark signal.

    If the dark signal does not depend linearly on exposure time, make
    sure the dark current CKD is valid for the used exposure time.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    dark_current
        Detector map of dark current [counts/s].

    """
    l1_product['proc_level'] = ProcLevel.noise
    l1_product['image'] += dark_current * l1_product['exptimes'][..., None]


def noise(l1_product: L1,
          ckd: CKDNoise,
          dark_current: npt.NDArray[np.float64],
          seed: int) -> None:
    """Add random noise to signal.

    The signal is changed.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Noise CKD consisting of maps of read_noise [counts] and the
        conversion gain [e/counts].
    dark_current
        Detector map of dark current [counts/s].
    seed
        Seed for random number generator.

    """
    l1_product['proc_level'] = ProcLevel.dark_offset
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = (ckd['read_noise']**2
                + l1_product['image'] / ckd['conversion_gain'])
    std = np.sqrt(np.clip(variance, 0, None))
    rng = np.random.default_rng(seed)
    l1_product['image'] += rng.normal(0.0, std, l1_product['image'].shape)


def dark_offset(l1_product: L1, offset: npt.NDArray[np.float64]) -> None:
    """Incorporate offset.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    offset
        Detector map of offset [counts].

    """
    l1_product['proc_level'] = ProcLevel.raw
    l1_product['image'] += offset


def coadding_and_binning(l1_product: L1,
                         binning_table: BinningTable) -> None:
    """Coadd over time and bin over the detector, both as sums.

    Coaddition effectively makes exact copies of existing frames and
    sums them, i.e. if noise was added, that is exactly the same in
    each frame copy.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    coad_factor
        Coaddition factor.
    binning_table
        Bin index of each pixel in an unbinned frame and number of
        pixels in each bin of a binned frame.

    """
    l1_product['proc_level'] = ProcLevel.l1a
    n_alt = l1_product['image'].shape[0]
    binned_images = np.zeros((n_alt, len(binning_table['count_table'])))
    for i_alt in range(l1_product['image'].shape[0]):
        for idx, idx_binned in enumerate(binning_table['bin_indices'].ravel()):
            binned_images[i_alt, idx_binned] += l1_product['image'][i_alt, idx]
    l1_product['image'] = binned_images
    l1_product['image'] = l1_product['image'] * l1_product['coad_factors'][0]
    # In reality, the noise realization is different in every read-out
    # and the signal is rounded before summing. Here the noise is not
    # calculated for each read-out separately. Then the best
    # approximation is to round after summing.
    l1_product['image'] = np.round(l1_product['image'])
