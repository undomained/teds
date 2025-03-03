# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Fixtures and settings used by more than one modules."""
from pathlib import Path
from pyquaternion import Quaternion
from scipy.fft import fft2
import pytest
import numpy as np

from teds.l1al1b.io import read_binning_table
from teds.l1al1b.types import CKDDark
from teds.l1al1b.types import CKDNoise
from teds.l1al1b.types import CKDNonlin
from teds.l1al1b.types import CKDPRNU
from teds.l1al1b.types import CKDRadiometric
from teds.l1al1b.types import CKDSpectral
from teds.l1al1b.types import CKDStray
from teds.l1al1b.types import CKDSwath
from teds.l1al1b.types import L1

_cur_dir = Path(__file__).parent
_fix_dir = _cur_dir / 'common'


# Fixtures used by the instrument model and L1A-L1B processor
@pytest.fixture
def binning_table(scop='session'):
    signal = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4')
    return read_binning_table('', 0, signal.shape[0], signal.shape[1])


@pytest.fixture
def signal_i4(scope='session'):
    signal = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4')
    return signal.reshape((1, -1))


@pytest.fixture
def signal(signal_i4, scope='session'):
    return signal_i4.astype('f8')


@pytest.fixture
def ckd_pixel_mask(scope='session'):
    return np.loadtxt(_fix_dir / 'ckd_pixel_mask.txt', '?').ravel()


@pytest.fixture
def ckd_dark(scope='session'):
    return CKDDark(np.loadtxt(_fix_dir / 'ckd_dark_offset.txt', 'f8').ravel(),
                   np.loadtxt(_fix_dir / 'ckd_dark_current.txt', 'f8').ravel())


@pytest.fixture
def ckd_noise(scope='session'):
    return CKDNoise(
        np.loadtxt(_fix_dir / 'ckd_conversion_gain.txt', 'f8').ravel(),
        np.loadtxt(_fix_dir / 'ckd_read_noise.txt', 'f8').ravel())


@pytest.fixture
def ckd_nonlin(scope='session'):
    data = np.loadtxt(_fix_dir / 'ckd_nonlin.txt', 'f8')
    return CKDNonlin(data[:, 0], data[:, 1])


@pytest.fixture
def ckd_prnu(scope='session'):
    return CKDPRNU(np.loadtxt(_fix_dir / 'ckd_prnu.txt', 'f8').ravel())


@pytest.fixture
def ckd_stray(scope='session'):
    # Dimensions and input data
    n_rows, n_cols = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4').shape
    kernel = np.loadtxt(_fix_dir / 'ckd_stray_kernel.txt', 'f8')
    # Get all variables for one kernel
    kernels_fft_1 = fft2(kernel)
    n_rows_fft, n_cols_fft = kernels_fft_1.shape
    eta = np.full((n_rows, n_cols), 0.1)
    weights_1 = np.ones(eta.shape)
    edges_1 = np.array([0, n_rows, 0, n_cols])
    # Duplicate for the second kernel
    kernels_fft = np.tile(kernels_fft_1, 2).reshape(2, n_rows_fft, n_cols_fft)
    weights = np.tile(weights_1, 2).reshape(2, n_rows, n_cols)
    edges = np.tile(edges_1, 2).reshape(2, -1)
    return CKDStray(kernels_fft, eta, weights, edges)


@pytest.fixture
def ckd_swath(scope='session'):
    return CKDSwath(
        np.loadtxt(_fix_dir / 'ckd_swath_act_angle.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_wavelength.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_act_map.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_wavelength_map.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_row_map.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_col_map.txt', 'f8'),
        np.loadtxt(_fix_dir / 'ckd_swath_los.txt', 'f8'))


@pytest.fixture
def ckd_spectral(scope='session'):
    return CKDSpectral(np.loadtxt(_fix_dir / 'ckd_spectral.txt', 'f8'))


@pytest.fixture
def ckd_radiometric(scope='session'):
    return CKDRadiometric(np.loadtxt(_fix_dir / 'ckd_radiometric.txt', 'f8'))


@pytest.fixture
def navigation_time(scope='session'):
    return np.loadtxt(_fix_dir / 'navigation_time.txt', 'f8')


@pytest.fixture
def navigation_orb_pos(scope='session'):
    return np.loadtxt(_fix_dir / 'navigation_orb_pos.txt', 'f8')


@pytest.fixture
def navigation_att_quat(scope='session'):
    return np.loadtxt(_fix_dir / 'navigation_att_quat.txt', 'f8')


@pytest.fixture
def navigation(navigation_time,
               navigation_orb_pos,
               navigation_att_quat,
               scope='session'):
    n_alt = len(navigation_time)
    att_quat = np.empty((n_alt,), dtype=Quaternion)
    for i in range(n_alt):
        q = navigation_att_quat[i, :]
        att_quat[i] = Quaternion(q[3], *q[:3]).normalised
    return Navigation(navigation_time,
                      navigation_orb_pos,
                      att_quat,
                      np.zeros(navigation_time.shape))


@pytest.fixture
def ckd(ckd_pixel_mask,
        ckd_dark,
        ckd_noise,
        ckd_nonlin,
        ckd_prnu,
        ckd_stray,
        ckd_swath,
        ckd_spectral,
        ckd_radiometric,
        scope='session'):
    n_rows, n_cols = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4').shape
    return CKD(n_rows,
               n_cols,
               ckd_pixel_mask,
               ckd_dark,
               ckd_noise,
               ckd_nonlin,
               ckd_prnu,
               ckd_stray,
               ckd_swath,
               ckd_spectral,
               ckd_radiometric)


@pytest.fixture
def l1(signal, scope='session'):
    l1 = L1.from_empty()
    l1.signal = signal
    l1.noise = np.ones(signal.shape)
    l1.spectra = np.loadtxt(_fix_dir / 'l1b_radiance.txt', 'f8')
    l1.spectra_noise = np.loadtxt(_fix_dir / 'l1b_radiance_stdev.txt', 'f8')
    l1.spectra = l1.spectra.reshape(
        (1, l1.spectra.shape[0], l1.spectra.shape[1]))
    l1.spectra_noise = l1.spectra_noise.reshape(
        (1, l1.spectra_noise.shape[0], l1.spectra_noise.shape[1]))
    l1.exposure_time = 0.0460833
    return l1
