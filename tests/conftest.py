# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Fixtures and settings used by more than one modules."""
from netCDF4 import Dataset
from pathlib import Path
from pyquaternion import Quaternion
from scipy.fft import fft2
from scipy.interpolate import CubicSpline
import pytest
import numpy as np

from teds.gm.types import Navigation
from teds.l1al1b.io import read_binning_table
from teds.l1al1b.binning import bin_data
from teds.l1al1b.types import BinningTable
from teds.l1al1b.types import CKD
from teds.l1al1b.types import CKDDark
from teds.l1al1b.types import CKDNoise
from teds.l1al1b.types import CKDNonlin
from teds.l1al1b.types import CKDPRNU
from teds.l1al1b.types import CKDRadiometric
from teds.l1al1b.types import CKDSpectral
from teds.l1al1b.types import CKDStray
from teds.l1al1b.types import CKDSwath
from teds.l1al1b.types import L1
from teds.sgm.atmosphere import Atmosphere

_cur_dir = Path(__file__).parent
_fix_dir = _cur_dir / 'common'


##################################################################
# Fixtures for instrument model (IM) and L1A-L1B processor (L1B) #
##################################################################

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
    kernels_fft = np.tile(kernels_fft_1.ravel(), 2).reshape(
        2, n_rows_fft, n_cols_fft)
    weights = np.tile(weights_1, 2).reshape(2, n_rows, n_cols)
    edges = np.tile(edges_1, 2).reshape(2, -1)
    return CKDStray(kernels_fft, eta, weights, edges)


@pytest.fixture
def ckd_swath(scope='session'):
    return CKDSwath(
        np.loadtxt(_fix_dir / 'ckd_swath_act_angle.txt', 'f8'),
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


@pytest.fixture
def sgm(l1, ckd_spectral, scope='session'):
    """Generate SGM line-by-line spectra."""
    n_alt, n_act, _ = l1.spectra.shape
    # Line-by-line (LBL) wavelength grid is based on the CKD
    # intermediate wavelength grid.
    lbl_wave_multiplier = 3
    margin = 1.0
    lbl_wavelengths = np.linspace(
        ckd_spectral.wavelengths[0] - margin,
        ckd_spectral.wavelengths[-1] + margin,
        lbl_wave_multiplier*len(ckd_spectral.wavelengths))
    lbl_spectra = np.empty((n_alt, n_act, len(lbl_wavelengths)))
    for i_alt in range(n_alt):
        for i_act in range(n_act):
            spline = CubicSpline(ckd_spectral.wavelengths,
                                 l1.spectra[i_alt, i_act, :])
            lbl_spectra[i_alt, i_act, :] = spline(lbl_wavelengths)
    l1.spectra = lbl_spectra
    l1.spectra_noise = np.ones(l1.spectra.shape)
    l1.wavelengths = lbl_wavelengths
    return l1


def gen_binning_table(n_rows, n_cols, bin_f):
    """Construct binning table based on detector dimensions and
    binning factor."""
    n_rows_red = n_rows // bin_f
    binned_indices = np.arange(
        0, n_rows_red * n_cols, 1, dtype=int).reshape(n_rows_red, n_cols)
    binning_table = np.zeros((n_rows, n_cols), dtype=int)
    for i_bin in range(bin_f):
        binning_table[i_bin:bin_f*n_rows_red:bin_f, :] = binned_indices
    bin_f_new = n_rows - bin_f * n_rows_red
    if bin_f_new > 0:
        next_idx = binned_indices.ravel()[-1] + 1
        for i_row in range(bin_f * n_rows_red, n_rows):
            binning_table[i_row, :] = np.arange(
                next_idx, next_idx + n_cols, dtype=int)
    uniq_bin_indices, tmp_count_table = np.unique(binning_table.ravel(),
                                                  return_counts=True)
    count_table = np.zeros(np.max(uniq_bin_indices) + 1, np.int32)
    count_table[uniq_bin_indices] = tmp_count_table
    return BinningTable(binning_table, count_table)


@pytest.fixture
def binningtable_file(tmp_path, scope='session'):
    """Fixture for writing binning table to temporary path.

    All binning tables are constructed on the fly and saved to the
    NetCDF file.

    Returns
    -------
        Full path to the temporary file.

    """
    filepath = tmp_path / 'tango_binningtable.nc'
    nc = Dataset(filepath, 'w')
    n_rows, n_cols = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4').shape
    nc.createDimension('row', n_rows)
    nc.createDimension('column', n_cols)
    for bin_f in range(1, 6):
        binning_table = gen_binning_table(n_rows, n_cols, bin_f)
        grp = nc.createGroup(f'Table_{bin_f}')
        grp.createDimension('bins', len(binning_table.count_table))
        var = grp.createVariable('binning_table', 'u4', ('row', 'column'))
        var[:] = binning_table.bin_indices
        var = grp.createVariable('count_table', 'u2', 'bins')
        var[:] = binning_table.count_table
    return filepath


@pytest.fixture
def l1a_file(request, tmp_path, navigation, l1, scope='session'):
    """Fixture for writing an L1A product to temporary path.

    Returns
    -------
        Full path to the temporary file.

    """
    filepath = tmp_path / 'tango_l1a.nc'
    nc = Dataset(filepath, 'w')

    # First bin L1A data
    binning_factor = request.param['bin'] if hasattr(request, 'param') else 1
    n_rows, n_cols = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4').shape
    binning_table = gen_binning_table(n_rows, n_cols, binning_factor)
    n_bins = len(binning_table.count_table)
    l1.binning_table_id = binning_factor
    signal_binned = np.zeros((l1.signal.shape[0], n_bins))
    for i_alt in range(l1.signal.shape[0]):
        signal_binned[i_alt, :] = bin_data(
            binning_table, l1.signal[i_alt, :], False)
    l1.signal = signal_binned

    nc.createDimension('along_track_sample', l1.signal.shape[0])
    nc.createDimension('bin', n_bins)
    nc.product_type = 'L1A'
    nc.instrument = 'Carbon'
    # Image attributes
    grp = nc.createGroup('image_attributes')
    var = grp.createVariable('time', 'f8', 'along_track_sample')
    var.units = 'seconds since 2025-01-01'
    var[:] = 44025.4149999619
    var = grp.createVariable('tai_seconds', 'u4', 'along_track_sample')
    var[:] = 2038997625
    var = grp.createVariable('tai_subsec', 'u2', 'along_track_sample')
    var[:] = 27197
    var = grp.createVariable('binning_table', 'i1')
    var[:] = l1.binning_table_id
    var = grp.createVariable('nr_coadditions', 'u2')
    var[:] = 2
    var = grp.createVariable('exposure_time', 'f8')
    var[:] = 0.01724385
    # Science data
    grp = nc.createGroup('science_data')
    var = grp.createVariable(
        'detector_image', 'i4', ('along_track_sample', 'bin'))
    var[:] = l1.signal
    # Navigation data
    grp = nc.createGroup('navigation_data')
    grp.createDimension('time', len(navigation.time))
    grp.createDimension('vector_elements', 3)
    grp.createDimension('quaternion_elements', 4)
    var = grp.createVariable('time', 'f8', 'time')
    var[:] = navigation.time
    var = grp.createVariable('orb_pos', 'f8', ('time', 'vector_elements'))
    var[:] = navigation.orb_pos
    var = grp.createVariable('att_quat', 'f8', ('time', 'quaternion_elements'))
    att_quat = np.empty((len(navigation.att_quat), 4))
    for i in range(len(navigation.att_quat)):
        att_quat[i, :] = np.roll(navigation.att_quat[i].elements, -1)
    var[:] = att_quat
    return filepath


@pytest.fixture
def ckd_file(tmp_path, ckd, scope='session'):
    """Fixture for writing CKD to temporary path.

    Returns
    -------
        Full path to the temporary file.

    """
    filepath = tmp_path / 'tango_ckd.nc'
    nc = Dataset(filepath, 'w')
    nc.createDimension('detector_row', ckd.n_detector_rows)
    nc.createDimension('detector_column', ckd.n_detector_cols)
    nc.createDimension('across_track_sample', len(ckd.swath.act_angles))
    nc.createDimension('wavelength', len(ckd.spectral.wavelengths))
    nc.createDimension('vector', 3)
    detector_shape = (ckd.n_detector_rows, ckd.n_detector_cols)
    dim_detector_shape = ('detector_row', 'detector_column')
    # Pixel mask
    var = nc.createVariable('pixel_mask', 'u1', dim_detector_shape)
    var[:] = ckd.pixel_mask.reshape(detector_shape)
    # Dark CKD
    grp = nc.createGroup('dark')
    var = grp.createVariable('offset', 'f8', dim_detector_shape)
    var[:] = ckd.dark.offset.reshape(detector_shape)
    var = grp.createVariable('current', 'f8', dim_detector_shape)
    var[:] = ckd.dark.current.reshape(detector_shape)
    # Noise CKD
    grp = nc.createGroup('noise')
    var = grp.createVariable('g', 'f8', dim_detector_shape)
    var[:] = ckd.noise.conversion_gain.reshape(detector_shape)
    var = grp.createVariable('n', 'f8', dim_detector_shape)
    var[:] = ckd.noise.read_noise.reshape(detector_shape)
    # Nonlinearity CKD
    grp = nc.createGroup('nonlinearity')
    grp.createDimension('knots', len(ckd.nonlin.expected))
    var = grp.createVariable('observed', 'f8', 'knots')
    var[:] = ckd.nonlin.observed
    var = grp.createVariable('expected', 'f8', 'knots')
    var[:] = ckd.nonlin.expected
    # PRNU CKD
    grp = nc.createGroup('prnu')
    var = grp.createVariable('prnu', 'f8', dim_detector_shape)
    var[:] = ckd.prnu.prnu_qe.reshape(detector_shape)
    # Stray light CKD
    grp = nc.createGroup('stray')
    n_kernels = ckd.stray.kernels_fft.shape[0]
    # Save kernels in non-packed format
    kernel_n_rows, kernel_n_cols = np.loadtxt(
        _fix_dir / 'ckd_stray_kernel.txt', 'f8').shape
    kernel_fft_size = kernel_n_rows * kernel_n_cols  # In units of 16 bytes
    grp.createDimension('kernel', n_kernels)
    # Total size is in units of 8 bytes
    grp.createDimension('kernels_fft_size', 2 * n_kernels * kernel_fft_size)
    grp.createDimension('edges_of_box', 4)
    var = grp.createVariable('kernel_rows', 'i4', 'kernel')
    var[:] = kernel_n_rows
    var = grp.createVariable('kernel_cols', 'i4', 'kernel')
    var[:] = kernel_n_cols
    var = grp.createVariable('kernel_fft_sizes', 'i4', 'kernel')
    var[:] = kernel_fft_size
    var = grp.createVariable('kernels_fft', 'f8', 'kernels_fft_size')
    for i_kern in range(n_kernels):
        n_comp = 2
        for i_comp in range(n_comp):
            if i_comp == 0:
                data = np.real(ckd.stray.kernels_fft[i_kern])
            else:
                data = np.imag(ckd.stray.kernels_fft[i_kern])
            i_start = (i_kern * n_comp + i_comp) * kernel_fft_size
            var[i_start:i_start+kernel_fft_size] = data
    var = grp.createVariable('eta', 'f8', dim_detector_shape)
    var[:] = ckd.stray.eta
    var = grp.createVariable(
        'weights', 'f8', ('kernel', 'detector_row', 'detector_column'))
    var[:] = ckd.stray.weights
    var = grp.createVariable('edges', 'i4', ('kernel', 'edges_of_box'))
    var[:] = ckd.stray.edges
    # Swath CKD
    grp = nc.createGroup('swath')
    var = grp.createVariable('act_angle', 'f8', 'across_track_sample')
    var[:] = ckd.swath.act_angles
    var = grp.createVariable('act_map', 'f8', dim_detector_shape)
    var[:] = ckd.swath.act_map
    var = grp.createVariable('wavelength_map', 'f8', dim_detector_shape)
    var[:] = ckd.swath.wavelength_map
    var = grp.createVariable(
        'row_map', 'f8', ('across_track_sample', 'wavelength'))
    var[:] = ckd.swath.row_map
    var = grp.createVariable(
        'col_map', 'f8', ('across_track_sample', 'wavelength'))
    var[:] = ckd.swath.col_map
    var = grp.createVariable(
        'line_of_sight', 'f8', ('across_track_sample', 'vector'))
    var[:] = ckd.swath.line_of_sights
    # Spectral CKD
    grp = nc.createGroup('spectral')
    var = grp.createVariable('wavelength', 'f8', 'wavelength')
    var[:] = ckd.spectral.wavelengths
    # Radiometric CKD
    grp = nc.createGroup('radiometric')
    var = grp.createVariable(
        'radiometric', 'f8', ('across_track_sample', 'wavelength'))
    var[:] = ckd.radiometric.rad_corr
    return filepath


@pytest.fixture
def sgm_file(tmp_path, sgm, scope='session'):
    """Fixture for writing CKD to temporary path.

    Returns
    -------
        Full path to the temporary file.

    """
    filepath = tmp_path / 'tango_sgm.nc'
    nc = Dataset(filepath, 'w')
    nc.product_type = 'SGM'
    nc.instrument = 'Carbon'
    nc.createDimension('along_track_sample', sgm.spectra.shape[0])
    nc.createDimension('across_track_sample', sgm.spectra.shape[1])
    nc.createDimension('wavelength', sgm.spectra.shape[2])
    var = nc.createVariable('wavelength', 'f8', 'wavelength')
    var[:] = sgm.wavelengths
    var = nc.createVariable(
        'radiance',
        'f8',
        ('along_track_sample', 'across_track_sample', 'wavelength'))
    var[:] = sgm.spectra
    return filepath


##############################################
# Fixtures for scene generation module (SGM) #
##############################################

@pytest.fixture
def atmosphere_path(scope='session'):
    return _fix_dir / 'prof.AFGL.US.std'


@pytest.fixture
def atmosphere(scope='session'):
    nlay = 20
    dzlay = 1000
    psurf = 101300
    nlev = nlay + 1
    zlay = (np.arange(nlay - 1, -1, -1) + 0.5) * dzlay
    zlev = np.arange(nlev - 1, -1, -1) * dzlay
    return Atmosphere(zlay, zlev, psurf, _fix_dir / 'prof.AFGL.US.std')
