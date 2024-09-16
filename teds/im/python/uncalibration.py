# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""List of processes applied to data by the instrument model.

Each process undoes or "uncalibrates" a step that is performed by the
L1A-L1B processor, gradually bringing the data level from L1B to L1A
or anywhere in between.

"""
from teds.l1al1b.python.binning import bin_data
from teds.l1al1b.python.calibration import convolve_with_all_kernels
from teds.l1al1b.python.io import read_binning_pattern
from teds.l1al1b.python.l1b_types import CKDNoise
from teds.l1al1b.python.l1b_types import CKDNonlin
from teds.l1al1b.python.l1b_types import CKDStray
from teds.l1al1b.python.l1b_types import ProcLevel
from teds.l1al1b.python.l1b_types import L1

from scipy.interpolate import CubicSpline
import math
import numpy as np
import numpy.typing as npt


def convolve_with_isrf(signals_in: npt.NDArray[np.float64],
                       w_in: npt.NDArray[np.float64],
                       grid_out: npt.NDArray[np.float64],
                       fwhm: float,
                       shape: float = 2.0) -> npt.NDArray[np.float64]:
    """Convolve spectra with an ISRF and reduce wavelength grid.

    Help function for change_wavelength_grid.

    Args:
      signals_in:
        Spectra with dimensions ALT, ACT, spectral sample.
      w_in:
        Wavelengths corresponding to input spectra (1D).
      grid_out:
        Wavelengths corresponding to output spectra with dimensions
        ACT, spectral sample.
      fwhm:
        Full width at half-maximum.
      shape:
        Shape parameter of the generalized normal distribution.  Value
        2 for Gauss (default), towards 1 for stronger wings and large
        values for a more blocky shape.

    Returns:
      Spectra convolved with ISRF with dimensions ALT, ACT, spectral
      sample.

    """
    signals = np.ascontiguousarray(np.transpose(signals_in, axes=(1, 0, 2)))
    n_act, n_alt, n_win = signals.shape
    n_wout = grid_out.shape[1]
    new_signals = np.empty((n_wout, n_act, n_alt))
    dw = np.empty(n_win)
    dw[0] = w_in[1]-w_in[0]
    dw[1:-1] = (w_in[2:]-w_in[:-2])/2
    dw[-1] = w_in[-1]-w_in[-2]
    dw = np.abs(dw)
    nrm = np.log(2)**(1/shape)/(fwhm*math.gamma(1+1/shape))
    for act in range(n_act):
        for i, w0 in enumerate(grid_out[act]):
            w_diff = np.abs(w_in-w0)
            relevant = w_diff <= 3*fwhm
            isrf_points = np.zeros(n_win)
            isrf_points[relevant] = nrm*2**(-(2*w_diff[relevant]/fwhm)**shape)
            isrf_points *= dw
            integral = np.sum(isrf_points)
            new_signals[i, act] = np.dot(signals[act], isrf_points)
            # if wavelength grid is fine enough, the integral over the kernel
            # is already 1 and normalization has no effect
            # not used in C++ code
            if integral > 0.1:
                new_signals[i, act] /= integral
    return np.transpose(new_signals, axes=(2, 1, 0))


def change_wavelength_grid(l1_products: L1,
                           grid_out: npt.NDArray[np.float64],
                           convolve: bool,
                           fwhm: float = 0.45,
                           shape: float = 2.0) -> None:
    """Change wavelength grid and optionally convolve with ISRF.

    The ISRF has a fixed shape as a function of wavelength and does not depend
    on pixel.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      grid_out:
        New wavelength grid [nm].
      convolve:
        Convolve with ISRF (True) or only interpolate (False)
      fwhm:
        Full width at half-maximum [nm] of ISRF.
      shape:
        ISRF shape parameter. Value 2 for Gauss (default), towards 1
        for stronger wings and large values for a more blocky shape.

    """
    w_in = l1_products['wavelengths']  # SGM data, so 1D
    if convolve:
        # convolve spectra with ISRF assuming w_in is monotonically
        # increasing or decreasing (checked while reading)
        # extrapolated values are close to zero assuming no bad values
        # in the spectra in C++ code fixed to Gauss
        l1_products['spectra'] = convolve_with_isrf(
            l1_products['spectra'], w_in, grid_out, fwhm, shape)
        # all values the same
        if np.allclose(l1_products['spectra_noise'],
                       l1_products['spectra_noise'][0, 0, 0],
                       equal_nan=True):
            l1_products['spectra_noise'] = (
                l1_products['spectra_noise'][:, :, :grid_out.shape[1]])
        else:  # should never be the case
            l1_products['spectra_noise'] = convolve_with_isrf(
                l1_products['spectra_noise'], w_in, grid_out, fwhm, shape)
    else:
        # Only interpolate spectra assuming w_in is monotonically
        # increasing or decreasing (checked while reading).
        # Extrapolated values are NaN. scipy.interpolate.interp1d is
        # deprecated.
        increasing = np.all(np.diff(w_in) > 0)
        if not increasing:
            w_in = np.flip(w_in, axis=-1)
            l1_products['spectra'] = np.flip(l1_products['spectra'], axis=-1)
            l1_products['spectra_noise'] = np.flip(
                l1_products['spectra_noise'], axis=-1)
        l1_products['spectra'] = np.array(
            [[np.interp(w_out, w_in, s_in, left=np.nan, right=np.nan)
              for w_out, s_in in zip(grid_out, signal)]
             for signal in l1_products['spectra']])
        l1_products['spectra_noise'] = np.array(
            [[np.interp(w_out, w_in, s_in, left=np.nan, right=np.nan)
              for w_out, s_in in zip(grid_out, signal)]
             for signal in l1_products['spectra_noise']])
    l1_products['wavelengths'] = grid_out
    l1_products['proc_level'] = ProcLevel.l1b


def convert_from_radiance(l1_products: L1,
                          rad_corr: npt.NDArray[np.float64],
                          exptime: float) -> None:
    """Convert from spectral photon radiance [nm-1 s-1 sr-1 m-2] to
    counts.

    The quantum efficiency [e/ph] is taken into account by the PRNU
    step.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      rad_corr:
        Radiance responsivity correction factor [nm-1 sr-1 m-2].
      exptime:
        Exposure time [s] if not already in l1_products.

    """
    # Should always be the case
    if np.isnan(l1_products['exptimes']).all():
        l1_products['exptimes'] = np.full_like(l1_products['exptimes'],
                                               exptime)
    if np.min(l1_products['exptimes']) <= 0:
        raise SystemExit("ERROR: exposure time cannot be "
                         f"{np.min(l1_products['exptimes'])} s")
    # mean_gain = np.mean(conversion_gain)
    mean_gain = 1  # in C++ code effectively 1
    l1_products['spectra'] *= (
        l1_products['exptimes'][:, None, None]/(mean_gain*rad_corr))
    # not yet in C++ code
    l1_products['spectra_noise'] *= (
        l1_products['exptimes'][:, None, None]/(mean_gain*rad_corr))
    l1_products['proc_level'] = ProcLevel.swath


def add_boundary_knots(spline: CubicSpline) -> None:
    """Add knots infinitesimally to the left and right of a cubic
    spline.

    Additional intervals are added to have zero 2nd and 3rd
    derivatives, and to maintain the first derivative from whatever
    boundary condition was selected. The spline is modified in place.

    """
    # determine the slope at the left edge
    leftx = spline.x[0]
    lefty = spline(leftx)
    leftslope = spline(leftx, nu=1)
    # add a new breakpoint just to the left and use the
    # known slope to construct the PPoly coefficients
    leftxnext = np.nextafter(leftx, leftx - 1)
    leftynext = lefty + leftslope*(leftxnext - leftx)
    leftcoeffs = np.array([0*leftynext, 0*leftynext, leftslope, leftynext])
    spline.extend(np.expand_dims(leftcoeffs, 1), np.r_[leftxnext])
    # repeat with additional knots to the right
    rightx = spline.x[-1]
    righty = spline(rightx)
    rightslope = spline(rightx, nu=1)
    rightxnext = np.nextafter(rightx, rightx + 1)
    rightynext = righty + rightslope * (rightxnext - rightx)
    rightcoeffs = np.array(
        [0*rightynext, 0*rightynext, rightslope, rightynext])
    spline.extend(np.expand_dims(rightcoeffs, 1), np.r_[rightxnext])


def scene_to_detector(scene_data: npt.NDArray[np.float64],
                      spectrum_rows: npt.NDArray[np.float64],
                      n_detrows: int) -> npt.NDArray[np.float64]:
    """Map data from a scene to the detector.

    Data from a scene are mapped to infinitely thin curves on the
    detector (instead of simply to rows due to spatial distortion or
    keystone) and then interpolated along columns (not along constant
    wavelengths) at integer row coordinates, using cubic interpolation
    and linear extrapolation.

    Args:
      scene_data:
        Data with dimensions ALT angle, ACT angle, wavelength.
      spectrum_rows:
        Central row index (float) of a given spectrum (i.e. of data at
        a given ACT angle) and at each column, 0 is halfway first row.
      n_detrows:
        Number of detector rows.

    Returns:
      Data with dimensions ALT angle, row, column.

    """
    # Assuming spectrum_rows is monotonically increasing or decreasing
    # per column, but all columns in the same direction (checked while
    # reading). NaNs are interpolated unless (almost) all data in a
    # column are NaNs. C++ code does extrapolation at last rows
    # incorrectly.
    n_alt = len(scene_data)
    n_detcols = spectrum_rows.shape[1]
    increasing = np.all(np.diff(spectrum_rows[:, 0]) > 0)
    det_data = np.empty((n_alt, n_detrows, n_detcols))
    for w in range(n_detcols):
        # Which row indices have good values in all frames (for the
        # given column)
        good = ~np.isnan(scene_data[:, :, w]).any(axis=0)
        if np.count_nonzero(good) < 2:
            det_data[:, :, w] = np.nan
        else:
            if increasing:
                cs = CubicSpline(spectrum_rows[good, w],
                                 scene_data[:, good, w],
                                 axis=1,
                                 bc_type='natural')
            else:
                cs = CubicSpline(spectrum_rows[good, w][::-1],
                                 scene_data[:, good, w][:, ::-1],
                                 axis=1,
                                 bc_type='natural')
            add_boundary_knots(cs)
            det_data[:, :, w] = cs(np.arange(n_detrows))
    return det_data


def map_to_detector(l1_products: L1,
                    spectrum_rows: npt.NDArray[np.float64],
                    n_detrows: int) -> None:
    """Map spectra to detector.

    Spectra are mapped to infinitely thin curves on the detector (instead of
    simply to rows due to spatial distortion or keystone) and then interpolated
    along columns (not along constant wavelengths) at integer row coordinates,
    using cubic interpolation and linear extrapolation.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      spectrum_rows:
        Central row index (float) of a given spectrum and at each
        column.
      n_detrows:
        Number of detector rows.

    """
    l1_products['signal'] = scene_to_detector(
        l1_products['spectra'],
        spectrum_rows,
        n_detrows).reshape(len(l1_products['spectra']), -1)
    # C++ code does not do l1_products['noise'] yet
    l1_products['noise'] = scene_to_detector(
        l1_products['spectra_noise'],
        spectrum_rows,
        n_detrows).reshape(len(l1_products['spectra_noise']), -1)
    l1_products['proc_level'] = ProcLevel.stray


def stray_light(l1_products: L1, ckd: CKDStray) -> None:
    """Add stray light to the image.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      ckd:
        Stray light CKD containing a list of the Fourier transforms of
        kernels, weights of subimages, and an 'edges' array which
        specifies the location of each subimage within the original
        image.

    """
    eta = ckd['eta'].reshape(l1_products['signal'][0, :].shape[0])
    for i_alt in range(l1_products['signal'].shape[0]):
        image = l1_products['signal'][i_alt]
        stray_convolved = convolve_with_all_kernels(image, ckd)
        image_convolved = (1 - eta) * image + stray_convolved
        l1_products['signal'][i_alt, :] = image_convolved


def include_prnu(l1_products: L1, prnu_qe: npt.NDArray[np.float64]) -> None:
    """Incorporate PRNU and quantum efficiency.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      prnu_qe:
        Detector map of PRNU times quantum efficiency (not the
        correction).

    """
    # Assuming bad pixels are NaN already
    l1_products['signal'] *= prnu_qe.ravel()
    l1_products['noise'] *= prnu_qe.ravel()  # not yet in C++ code
    l1_products['proc_level'] = ProcLevel.nonlin


def include_nonlinearity(l1_products: L1, ckd: CKDNonlin) -> None:
    """Incorporate non-linearity.

    Only a dependence on signal is taken into account, not on
    intensity, so make sure the non-linearity CKD is valid for the
    used exposure time.

    The model is a linear interpolation of a set of points, so make
    sure the number of points is large enough to produce a smooth
    curve. Input data outside the model range are clipped. The model
    does not depend on pixel.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      ckd:
        Nonlinearity CKD consisting of an expected, linear signal
        [counts], and an observed, non-linear signal [counts].

    """
    # noise is not yet in C++ code
    dx = 0.001*np.min(np.diff(ckd['expected']))
    l1_products['noise'] *= (
        (np.interp(l1_products['signal']+dx,
                   ckd['expected'],
                   ckd['observed'])
         - np.interp(l1_products['signal']-dx,
                     ckd['expected'],
                     ckd['observed']))
        / (2*dx))
    l1_products['signal'] = np.interp(
        l1_products['signal'], ckd['expected'], ckd['observed'])
    l1_products['proc_level'] = ProcLevel.noise


def include_darksignal(l1_products: L1,
                       dark_current: npt.NDArray[np.float64]) -> None:
    """Incorporate dark signal.

    If the dark signal does not depend linearly on exposure time, make
    sure the dark current CKD is valid for the used exposure time.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      dark_current:
        Detector map of dark current [counts/s].

    """
    # Assuming bad pixels are NaN already
    l1_products['signal'] += (dark_current.ravel()
                              * l1_products['exptimes'][..., None])
    # Uncertainty of dark signal is not implemented, no change of
    # processing level.


def include_noise(l1_products: L1,
                  ckd: CKDNoise,
                  dark_current: npt.NDArray[np.float64],
                  seed: int) -> None:
    """Add random noise to signal.

    The signal is changed.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      ckd:
        Noise CKD consisting of maps of read_noise [counts] and the
        conversion gain [e/counts].
      dark_current:
        Detector map of dark current [counts/s].
      seed:
        Seed for random number generator.

    """
    # Assuming bad pixels are NaN already. Same seed as in C++ code
    # still gives different values due to different random generator.
    # Potential problem: if processing to 'dark' towards L1A, dark
    # signal has not been included yet. C++ code does not clip before
    # sqrt and does not take into account bad pixels the signal is
    # split in contributions from light and dark in case the dark
    # signal is negative.
    dark_signal = dark_current.ravel()*l1_products['exptimes'][..., None]
    photoresponse_signal = l1_products['signal'] - dark_signal
    # The absolute value of dark_signal should be taken because a
    # negative signal still increases the noise.
    variance = (ckd['read_noise'].ravel()**2
                + (photoresponse_signal + dark_signal)
                / ckd['conversion_gain'].ravel())
    std = np.sqrt(np.clip(variance, 0, None))
    rng = np.random.default_rng(seed)
    l1_products['signal'] += rng.normal(0.0, std, l1_products['signal'].shape)
    l1_products['proc_level'] = ProcLevel.dark


def include_offset(l1_products: L1, offset: npt.NDArray[np.float64]) -> None:
    """Incorporate offset.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      offset:
        Detector map of offset [counts].

    """
    # assuming bad pixels are NaN already
    l1_products['signal'] += offset.ravel()
    # binning factors for consistent RAW data, not in C++ code
    l1_products['noise'] = np.ones_like(l1_products['signal'])
    l1_products['proc_level'] = ProcLevel.raw


def include_coadding_and_binning(l1_products: L1,
                                 new_coad_factor: int,
                                 binning_file: str,
                                 binning_table_id: int,
                                 n_detrows: int,
                                 n_detcols: int) -> None:
    """Coadd over time and bin over the detector, both as sums.

    Coaddition effectively makes exact copies of existing frames and
    sums them, i.e. if noise was added, that is exactly the same in
    each frame copy.

    Args:
      l1_products:
        L1 products (signal and detector settings).
      new_coad_factor:
        Coaddition factor.
      binning_file:
        Path of file with binning patterns.
      binning_table_id:
        Identifier of binning pattern, 0 is unbinned.
      n_detrows:
        Number of detector rows.
      n_detcols:
        Number of detector columns.

    """
    # Assuming bad pixels are NaN already. C++ code does not take into
    # account bad pixels, so there CKD pixel mask has to be binned to
    # know which values are bad.
    l1_products['signal'] = l1_products['signal']*new_coad_factor
    # In C++ code performed at any step if not towards_l1b
    l1_products['coad_factors'] = np.full_like(
        l1_products['coad_factors'], new_coad_factor)
    new_bin_indices, new_count_table = read_binning_pattern(
        binning_file, binning_table_id, n_detrows, n_detcols)
    if binning_table_id > 0:
        l1_products['signal'] = bin_data(
            l1_products['signal'], new_bin_indices, new_count_table)
        # In C++ code performed at any step if not towards_l1b
        l1_products['binning_table_ids'] = np.full_like(
            l1_products['binning_table_ids'], binning_table_id)
    # C++ code does not clip
    l1_products['signal'] = np.clip(l1_products['signal'], 0, 2**31-1)
    # In reality, the noise realization is different in every read-out
    # and the signal is rounded before summing. Here the noise is not
    # calculated for each read-out separately. Then the best
    # approximation is to round after summing.
    l1_products['signal'] = np.round(l1_products['signal'])
    l1_products['proc_level'] = ProcLevel.l1a
