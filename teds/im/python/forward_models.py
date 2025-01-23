# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of processes applied to data by the instrument model.

Each process undoes or "uncalibrates" a step that is performed by the
L1A-L1B processor, gradually bringing the data level from L1B to L1A
or anywhere in between.

"""
from netCDF4 import Dataset
from scipy.interpolate import CubicSpline
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
    wavelengths_out
        New wavelength grid [nm].
    convolve
        Convolve with ISRF (True) or only interpolate (False)
    fwhm
        Full width at half-maximum [nm] of ISRF.
    shape
        ISRF shape parameter. Value 2 for Gauss (default), towards 1
        for stronger wings and large values for a more blocky shape.
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
    conv = np.empty((n_alt, n_act, len(wavelengths_out)))
    if convolve:
        # Convolve spectra with ISRF assuming wavelengths_in is
        # monotonically increasing or decreasing (checked while
        # reading). Extrapolated values are close to zero assuming no
        # bad values in the spectra.
        kernel = convolution.Kernel(
            l1_product.wavelengths, wavelengths_out, fwhm, shape)
        # If the spectra are not in memory then read them one by one
        # for the convolution.
        in_memory = l1_product.spectra.shape[-1] > 0
        if not in_memory:
            nc = Dataset(sgm_filename)
        for i_alt in tqdm(range(n_alt), unit=' ALT'):
            for i_act in range(n_act):
                if in_memory:
                    conv[i_alt, i_act, :] = kernel.convolve(
                        l1_product.spectra[i_alt, i_act, :])
                else:
                    conv[i_alt, i_act, :] = kernel.convolve(
                        nc['radiance'][alt_beg + i_alt, i_act, :].data)
        l1_product.spectra = conv
    else:
        for i_alt in tqdm(range(n_alt), unit=' ALT'):
            for i_act in range(n_act):
                conv[i_alt, i_act, :] = np.interp(
                    wavelengths_out,
                    l1_product.wavelengths,
                    l1_product.spectra[i_alt, i_act, :])
    l1_product.spectra = conv
    l1_product.wavelengths = np.tile(wavelengths_out, n_act).reshape(n_act, -1)


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
    l1_product.proc_level = ProcLevel.swath
    for i_alt in tqdm(range(l1_product.spectra.shape[0])):
        l1_product.spectra[i_alt, :, :] *= (
            l1_product.exposure_time / (rad_corr[0, 0]))


def map_to_detector(l1_product: L1,
                    ckd: CKDSwath,
                    wavelengths: npt.NDArray[np.float64],
                    exact_drawing: bool) -> None:
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
    exact_drawing
        Option to draw each spectrum to the detector using the "up"
        and "down" pixels. Use this for experimenting but not in
        production.

    """
    l1_product.proc_level = ProcLevel.stray
    n_alt = l1_product.spectra.shape[0]
    n_rows, n_cols = ckd.act_map.shape
    l1_product.signal = np.empty((n_alt, n_rows * n_cols))
    act_wavelength_map = np.column_stack((ckd.act_map.ravel(),
                                          ckd.wavelength_map.ravel()))
    for i_alt in tqdm(range(n_alt)):
        l1_product.signal[i_alt, :] = interpn(
            (ckd.act_angles, ckd.wavelengths),
            l1_product.spectra[i_alt, :, :],
            act_wavelength_map,
            method='quintic',
            bounds_error=False,
            fill_value=None).reshape(n_rows * n_cols)
    if not exact_drawing:
        return
    if l1_product.binning_table_id > 1:
        raise SystemExit(
            'error: exact drawing algorithm only works with binning 1x1')
    l1_product.signal[:] = 0
    # When using the exact drawing algorithm, first regrid row_map and
    # spectra from intermediate wavelengths to wavelengths
    # corresponding to detector columns.
    n_act = ckd.row_map.shape[0]
    row_map = np.empty((n_act, n_cols))
    for i_act in range(n_act):
        row_map[i_act, :] = CubicSpline(
            ckd.wavelengths, ckd.row_map[i_act, :])(wavelengths[i_act, :])
    for i_alt in tqdm(range(n_alt)):
        spectra = np.empty((n_act, n_cols))
        signal = l1_product.signal[i_alt, :].reshape((n_rows, n_cols))
        for i_act in range(n_act):
            spectra[i_act, :] = CubicSpline(
                ckd.wavelengths,
                l1_product.spectra[i_alt, i_act, :])(wavelengths[i_act, :])
        for i_act in range(n_act):
            for i_wave in range(n_cols):
                row_dn = int(row_map[i_act, i_wave])
                row_up = row_dn + 1
                weight = row_map[i_act, i_wave] - row_dn
                if abs(signal[row_dn, i_wave]) < 1e-100:
                    signal[row_dn, i_wave] = spectra[i_act, i_wave]
                signal[row_up, i_wave] = (
                    spectra[i_act, i_wave]
                    - weight * signal[row_dn, i_wave]) / (1 - weight)


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
    l1_product.proc_level = ProcLevel.noise
    l1_product.signal += dark_current * l1_product.exposure_time


def noise(l1_product: L1,
          ckd: CKDNoise,
          dark_current: npt.NDArray[np.float64],
          artificial_scaling: float,
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
    variance = (ckd.read_noise**2 + l1_product.signal / ckd.conversion_gain)
    std = artificial_scaling * np.sqrt(np.clip(variance, 0, None))
    rng = np.random.default_rng(seed)
    l1_product.signal += rng.normal(0.0, std, l1_product.signal.shape)


def dark_offset(l1_product: L1, offset: npt.NDArray[np.float64]) -> None:
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
    l1_product.proc_level = ProcLevel.l1a
    n_alt = l1_product.signal.shape[0]
    binned_signals = np.zeros((n_alt, len(binning_table.count_table)))
    for i_alt in range(l1_product.signal.shape[0]):
        for idx, idx_binned in enumerate(binning_table.bin_indices.ravel()):
            binned_signals[i_alt, idx_binned] += l1_product.signal[i_alt, idx]
    l1_product.signal = binned_signals
    l1_product.signal = l1_product.signal * l1_product.coad_factor
    # In reality, the noise realization is different in every read-out
    # and the signal is rounded before summing. Here the noise is not
    # calculated for each read-out separately. Then the best
    # approximation is to round after summing.
    l1_product.signal = np.round(l1_product.signal)
