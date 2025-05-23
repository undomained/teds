# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of processes applied to data by the instrument model.

Each process undoes or "uncalibrates" a step that is performed by the
L1A-L1B processor, gradually bringing the data level from L1B to L1A
or anywhere in between.

"""
from netCDF4 import Dataset
from scipy.interpolate import interpn
from tqdm import tqdm
import numpy as np
import numpy.typing as npt

from .types import Stokes
from teds.l1al1b.calibration import convolve_with_all_kernels
from teds.l1al1b.types import BinningTable
from teds.l1al1b.types import CKDNoise
from teds.l1al1b.types import CKDNonlin
from teds.l1al1b.types import CKDStray
from teds.l1al1b.types import CKDSwath
from teds.l1al1b.types import L1
from teds.l1al1b.types import ProcLevel
from teds.lib.convolution import Kernel


def apply_mueller(l1_product: L1,
                  stokes: Stokes,
                  mueller: npt.NDArray[np.floating]) -> None:
    """Add Mueller matrix elements to Stokes I,Q,U parameters.

    Compute I_pol = I + M_01 / M_00 * Q + M_02 / M_00 * U.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings)
    stokes
        Stokes Q and U parameters
    mueller
        Mueller matrix elements for each ALT, ACT, and wavelength

    """
    M_00 = mueller[0]
    M_01 = mueller[1]
    M_02 = mueller[2]
    l1_product.spectra += M_01 / M_00 * stokes.Q + M_02 / M_00 * stokes.U


def apply_isrf(l1_product: L1,
               isrf: Kernel,
               convolve: bool,
               sgm_filename: str,
               alt_beg: int) -> None:
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
    isrf
        Convolution kernel
    convolve
        Convolve with ISRF (True) or only interpolate (False)
    sgm_filename
        SGM radiance filename. Only required if the spectra are not in
        memory. Then read them one by one for ISRF convolution.
    alt_begin
        First ALT position for convolution if spectra are read from
        the SGM file. The last position can be derived from the
        dimensions of the allocated spectra variable.

    """
    l1_product.proc_level = ProcLevel.l1b
    n_alt = l1_product.spectra.shape[0]
    n_act = l1_product.spectra.shape[1]
    conv = np.empty((n_alt, n_act, len(isrf.wavelengths_out)))
    if convolve:
        # If the spectra are not in memory then read them one by one
        # for the convolution.
        in_memory = l1_product.spectra.shape[-1] > 0
        if not in_memory:
            nc = Dataset(sgm_filename)
        for i_alt in tqdm(range(n_alt), unit=' ALT'):
            for i_act in range(n_act):
                if in_memory:
                    conv[i_alt, i_act, :] = isrf.convolve(
                        l1_product.spectra[i_alt, i_act, :])
                else:
                    conv[i_alt, i_act, :] = isrf.convolve(
                        nc['radiance'][alt_beg + i_alt, i_act, :].data)
        l1_product.spectra = conv
    else:
        for i_alt in tqdm(range(n_alt), unit=' ALT'):
            for i_act in range(n_act):
                conv[i_alt, i_act, :] = np.interp(
                    isrf.wavelengths_out,
                    l1_product.wavelengths,
                    l1_product.spectra[i_alt, i_act, :])
    l1_product.spectra = conv
    l1_product.wavelengths = isrf.wavelengths_out


def radiometric(l1_product: L1, rad_corr: npt.NDArray[np.floating]) -> None:
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
    l1_product.proc_level = ProcLevel.swath
    for i_alt in tqdm(range(l1_product.spectra.shape[0])):
        l1_product.spectra[i_alt, :, :] *= (
            l1_product.exposure_time / (rad_corr[0, 0]))


def map_to_detector(l1_product: L1,
                    ckd: CKDSwath,
                    wavelengths: npt.NDArray[np.floating],
                    b_spline_order: int) -> None:
    """Map spectra to detector.

    Spectra are mapped to infinitely thin curves on the detector (instead of
    simply to rows due to spatial distortion or keystone) and then interpolated
    along columns (not along constant wavelengths) at integer row coordinates,
    using cubic interpolation and linear extrapolation.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Swath section of the CKD
    wavelengths
        Main CKD wavelength grids (not the intermediate grids). Only
        used if exact_drawing is true.
    b_spline_order
        2D B-spline order

    """
    l1_product.proc_level = ProcLevel.stray
    n_alt = l1_product.spectra.shape[0]
    n_rows, n_cols = ckd.act_map.shape
    l1_product.signal = np.empty((n_alt, n_rows * n_cols))
    act_wavelength_map = np.column_stack((ckd.act_map.ravel(),
                                          ckd.wavelength_map.ravel()))
    b_spline_methods = {1: 'linear', 3: 'cubic', 5: 'quintic'}
    for i_alt in tqdm(range(n_alt)):
        l1_product.signal[i_alt, :] = interpn(
            (ckd.act_angles, wavelengths),
            l1_product.spectra[i_alt, :, :],
            act_wavelength_map,
            method=b_spline_methods[b_spline_order],
            bounds_error=False,
            fill_value=None).reshape(n_rows * n_cols)


def stray_light(l1_product: L1, ckd: CKDStray) -> None:
    """Add stray light to the signal.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    ckd
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subsignals, and an 'edges' array which
        specifies the location of each subsignal within the original
        signal.

    """
    l1_product.proc_level = ProcLevel.prnu
    n_alt = l1_product.signal.shape[0]
    eta = ckd.eta.ravel()
    for i_alt in tqdm(range(n_alt)):
        signal = l1_product.signal[i_alt]
        stray_convolved = convolve_with_all_kernels(signal, ckd)
        signal_convolved = (1 - eta) * signal + stray_convolved
        l1_product.signal[i_alt, :] = signal_convolved


def prnu(l1_product: L1, prnu_qe: npt.NDArray[np.floating]) -> None:
    """Incorporate PRNU and quantum efficiency.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    prnu_qe
        Detector map of PRNU times quantum efficiency (not the
        correction).

    """
    l1_product.proc_level = ProcLevel.nonlin
    l1_product.signal *= prnu_qe


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
    l1_product.proc_level = ProcLevel.dark_current
    l1_product.signal = np.interp(
        l1_product.signal, ckd.expected, ckd.observed)


def dark_current(l1_product: L1,
                 dark_current: npt.NDArray[np.floating]) -> None:
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
    l1_product.proc_level = ProcLevel.noise
    l1_product.signal += dark_current * l1_product.exposure_time


def noise(l1_product: L1,
          ckd: CKDNoise,
          dark_current: npt.NDArray[np.floating],
          n_coadditions: int,
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
    l1_product.proc_level = ProcLevel.dark_offset
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = ((ckd.read_noise**2 + l1_product.signal / ckd.conversion_gain)
                / n_coadditions)
    std = np.sqrt(np.clip(variance, 0, None))
    rng = np.random.default_rng(seed)
    l1_product.signal += rng.normal(0.0, std, l1_product.signal.shape)


def dark_offset(l1_product: L1, offset: npt.NDArray[np.floating]) -> None:
    """Incorporate offset.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    offset
        Detector map of offset [counts].

    """
    l1_product.proc_level = ProcLevel.raw
    l1_product.signal += offset


def bin_detector_images(l1_product: L1,
                        binning_table_id: int,
                        binning_table: BinningTable) -> None:
    """Bin all detector images.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    binning_table_id
        Binning factor
    binning_table
        Bin index of each pixel in an unbinned frame and number of
        pixels in each bin of a binned frame.

    """
    n_alt = l1_product.signal.shape[0]
    binned_signals = np.zeros((n_alt, len(binning_table.count_table)))
    for i_alt in range(l1_product.signal.shape[0]):
        for idx, idx_binned in enumerate(binning_table.bin_indices.ravel()):
            binned_signals[i_alt, idx_binned] += l1_product.signal[i_alt, idx]
    l1_product.signal = binned_signals
    l1_product.binning_table_id = binning_table_id


def coadd_and_adc(l1_product: L1, n_coadditions: int) -> None:
    """Coadd and convert detector images to integer format.

    Coaddition effectively makes exact copies of existing frames and
    sums them, i.e. if noise was added, that is exactly the same in
    each frame copy.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    n_coadditions
        Coaddition factor.

    """
    l1_product.proc_level = ProcLevel.l1a
    l1_product.coad_factor = n_coadditions
    l1_product.signal = l1_product.signal * l1_product.coad_factor
    # In reality, the noise realization is different in every read-out
    # and the signal is rounded before summing. Here the noise is not
    # calculated for each read-out separately. Then the best
    # approximation is to round after summing.
    l1_product.signal = np.round(l1_product.signal)
