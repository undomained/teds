# This source code is licensed under the 3-clause BSD license found
# in the LICENSE file in the root directory of this project.
"""CKD generator for Tango Carbon

This module compiles a calibration key data (CKD) file using
experimental or simulated calibration data as input. Depending on the
CKD step, the data could be copied directly from an input file (raw
CKD) or undergo some modification first such as interpolation to a
different temperature or unit conversion.

Usage:

>>> from teds.gm import geometry_module
>>> geometry_module(yaml.safe_load(open('gm.yaml')))

"""
from netCDF4 import Dataset
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import RBFInterpolator
from scipy.interpolate import interpn
import pandas as pd
import numpy as np
import numpy.typing as npt

from teds import log
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info
import teds.l1al1b.python.types as tp


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration parameters.

    """
    ckd_file = Path(config['io_files']['ckd'])
    if not ckd_file.is_file():
        raise SystemExit(f"ERROR: {ckd_file} not found")


def gen_dark(conf: dict) -> tp.CKDDark:
    """Generate detector dark CKD.

    Parameters
    ----------
    conf
        settings read from the YAML file

    Returns
    -------
        CKD containing dark offset and current.

    """
    nc = Dataset(conf['io_files']['detector'])
    spline = PchipInterpolator(
        nc['temperature'][:], np.nan_to_num(nc['offset'][:]), 0)
    offset = spline(conf['temperature'])
    spline = PchipInterpolator(
        nc['temperature'][:], np.nan_to_num(nc['dark_current'][:]), 0)
    current = spline(conf['temperature'])
    return tp.CKDDark(offset, current)


def gen_noise(conf: dict) -> tp.CKDNoise:
    """Generate detector noise CKD."""
    nc = Dataset(conf['io_files']['detector'])
    spline = CubicSpline(nc['temperature'][:], 1 / nc['conversion_gain'][:])
    gain = spline(conf['temperature'])
    spline = PchipInterpolator(
        nc['temperature'][:], np.nan_to_num(nc['read_noise'][:]), 0)
    read_noise = spline(conf['temperature'])
    return tp.CKDNoise(gain, read_noise)


def gen_nonlin(conf: dict) -> tp.CKDNonlin:
    """Generate detector nonlinearity CKD."""
    # Dummy implementation for now
    n_knots = 100
    observed = np.linspace(-1000, 10000, n_knots)
    expected = np.linspace(-1000, 10000, n_knots)
    return tp.CKDNonlin(observed, expected)


def gen_spectral(conf: dict) -> tp.CKDSpectral:
    """Generate L1B wavelength grid."""
    wavelengths = np.linspace(
        conf['wavelength_min'], conf['wavelength_max'], conf['n_wavelengths'])
    return tp.CKDSpectral(wavelengths)


def gen_swath(conf: dict, ckd_spectral: tp.CKDSpectral) -> tp.CKDSwath:
    """Generate swath CKD."""
    l1b_act_angles = np.linspace(conf['l1b_act_angle_min'],
                                 conf['l1b_act_angle_max'],
                                 conf['n_l1b_act_angles'])

    # STEP 1 - Spot measurement data.

    # Spot across track angles, deg
    spot_act = np.asarray(pd.read_csv(
        open(conf['spot_act_angle_file']), sep='\t', header=None))[:, 0]
    # Spot wavelengths, nm
    spot_wavelengths = np.asarray(pd.read_csv(
        open(conf['spot_wavelength_file']), sep='\t', header=None))[:, 0]
    # Spot distances from the center detector row, mm
    row_distances = pd.read_csv(
        open(conf['spot_row_distance_file']), sep='\t', header=None)
    # Spot distances from the center detector column, mm
    col_distances = pd.read_csv(
        open(conf['spot_col_distance_file']), sep='\t', header=None)
    # Convert distances to fractional row and column indices of the pixels
    pixel_size = 0.015  # mm
    row_indices = np.asarray(row_distances) / pixel_size
    col_indices = np.asarray(col_distances) / pixel_size
    if not conf['enable_keystone']:
        log.info('  Switching off keystone')
        for i_col in range(row_indices.shape[1]):
            row_indices[:, i_col] = row_indices[:, 0]
    if not conf['enable_smile']:
        log.info('  Switching off spectral smile')
        for i_row in range(col_indices.shape[0]):
            col_indices[i_row, :] = col_indices[0, :]
    # Spot measurement data assumes a 512x640 detector. For any other
    # dimensions the generated pixel indices need to be scaled
    # appropriately.
    row_indices *= conf['n_detector_rows'] / 512
    col_indices *= conf['n_detector_cols'] / 640
    # Shift to the detector center
    row_indices += conf['n_detector_rows'] / 2
    col_indices += conf['n_detector_cols'] / 2

    # STEP 2 - ACT & wavelength mapping to pixels

    # All row and column indices used in spot measurements
    spot_rows_cols = np.column_stack((row_indices.ravel(),
                                      col_indices.ravel()))
    # ACT angles of all spot measurements. Like spot_act but values
    # are repeated across wavelengths.
    spot_act_all = np.repeat(spot_act, len(spot_wavelengths))
    # Wavelengths of all spot measurements. Like spot_wavelengths but
    # values are repeated across ACT angles.
    spot_wavelengths_all = np.tile(spot_wavelengths, len(spot_act))
    # Regular detector grid
    det_grid = np.mgrid[:conf['n_detector_rows'],
                        :conf['n_detector_cols']].reshape(2, -1).T
    # Interpolate spot measurement ACT angles to all pixels (to all
    # integer row and column indices).
    act_map = RBFInterpolator(
        spot_rows_cols,
        spot_act_all,
        kernel='cubic')(det_grid).reshape(conf['n_detector_rows'],
                                          conf['n_detector_cols'])
    # Similarly, interpolate spot measurement wavelengths to all pixels
    wavelength_map = RBFInterpolator(
        spot_rows_cols,
        spot_wavelengths_all,
        kernel='cubic')(det_grid).reshape(conf['n_detector_rows'],
                                          conf['n_detector_cols'])

    # STEP 3 - Row & column index mapping to L1B ACTs & wavelengths

    # ACT angles and wavelengths of all pixels. In L1B, we need to
    # interpolate from this to the target ACT and wavelength values.
    act_wavelength_map = np.column_stack((act_map.ravel(),
                                          wavelength_map.ravel()))
    # Row of each detector pixel
    det_rows = np.repeat(np.arange(conf['n_detector_rows']),
                         conf['n_detector_cols'])
    # Column of each detector pixel
    det_cols = np.tile(np.arange(conf['n_detector_cols']),
                       conf['n_detector_rows'])
    # Grid of target L1B act angles and wavelengths
    l1b_grid = np.array(np.meshgrid(l1b_act_angles,
                                    ckd_spectral.wavelengths,
                                    indexing='ij')).reshape(2, -1).T
    # Row index of each L1B spectrum
    log.info('  Generating row index of each L1B spectrum')
    row_map = RBFInterpolator(
        act_wavelength_map,
        det_rows,
        kernel='cubic',
        neighbors=49)(l1b_grid).reshape(len(l1b_act_angles),
                                        len(ckd_spectral.wavelengths))

    # Column index of each L1B spectrum
    log.info('  Generating column index of each L1B spectrum')
    col_map = RBFInterpolator(
        act_wavelength_map,
        det_cols,
        kernel='cubic',
        neighbors=49)(l1b_grid).reshape(len(l1b_act_angles),
                                        len(ckd_spectral.wavelengths))

    # Line-of-sight vectors
    l1b_act_angles_rad = np.deg2rad(l1b_act_angles)
    los = np.stack((np.zeros(len(l1b_act_angles_rad)),
                    np.sin(l1b_act_angles_rad),
                    np.cos(l1b_act_angles_rad)), 1)

    return tp.CKDSwath(l1b_act_angles,
                       act_map,
                       wavelength_map,
                       row_map,
                       col_map,
                       los)


def gen_prnu(conf: dict,
             ckd_swath: tp.CKDSwath,
             ckd_spectral: tp.CKDSpectral) -> tp.CKDPRNU:
    """Generate PRNU x QE CKD."""
    nc = Dataset(conf['io_files']['detector'])
    # QE is more difficult because it is provided per wavelength.
    # First step is to interpolate QE onto the target temperature.
    qe = np.zeros(nc.dimensions['wavelength'].size)
    temperatures_qe = nc['temperature_qe'][:]
    for i in range(len(qe)):
        spline = CubicSpline(temperatures_qe,
                             nc['quantum_efficiency'][:, i])
        qe[i] = spline(conf['temperature'])
    # Find the shortest CKD wavelength and pad the QE values so there
    # is always data at all wavelengths.
    min_ckd_wavelength = ckd_spectral.wavelengths.min() - 1.0  # for safety
    qe_wavelengths = nc['wavelength'][:].data
    qe_wavelength_step = qe_wavelengths[1] - qe_wavelengths[0]
    while qe_wavelengths[0] > min_ckd_wavelength:
        next_wavelength = qe_wavelengths[0] - qe_wavelength_step
        next_qe = qe[0]
        qe_wavelengths = np.insert(qe_wavelengths, 0, next_wavelength)
        qe = np.insert(qe, 0, next_qe)
    # Next interpolate QE onto the L1B spectra
    qe_spline = CubicSpline(qe_wavelengths, qe)
    qe_spectra = np.zeros((len(ckd_swath.act_angles),
                           len(ckd_spectral.wavelengths)))
    qe_spectrum = qe_spline(ckd_spectral.wavelengths)
    qe_spectra = np.tile(qe_spectrum,
                         qe_spectra.shape[0]).reshape(qe_spectra.shape)
    # At this point, we have a QE value for each L1B spectrum
    act_wave_map_out = np.column_stack((ckd_swath.act_map.ravel(),
                                        ckd_swath.wavelength_map.ravel()))
    qe_map = interpn(
            (ckd_swath.act_angles, ckd_spectral.wavelengths),
            qe_spectra,
            act_wave_map_out,
            method='cubic',
            bounds_error=False,
            fill_value=None).reshape(ckd_swath.act_map.shape)
    prnu_in = np.nan_to_num(nc['prnu'][:].astype('f8'))
    spline = PchipInterpolator(nc['temperature'][:], prnu_in, 0)
    prnu = spline(conf['temperature'])
    # Extrapolate over bad values near the detector edges
    margin = 2
    row_beg = margin
    row_end = prnu.shape[0] - margin
    col_beg = margin
    col_end = prnu.shape[1] - margin
    prnu = prnu[row_beg:row_end, col_beg:col_end]
    det_grid = np.mgrid[
        :conf['n_detector_rows'], :conf['n_detector_cols']].reshape(2, -1).T
    prnu = interpn((np.arange(row_beg, row_end), np.arange(col_beg, col_end)),
                   prnu,
                   det_grid,
                   method='linear',
                   bounds_error=False,
                   fill_value=None).reshape((conf['n_detector_rows'],
                                             conf['n_detector_cols']))
    # Combine PRNU and QE into one variable (for now)
    prnu = prnu * qe_map
    return tp.CKDPRNU(prnu)


def gen_radiometric(conf: dict,
                    ckd_swath: tp.CKDSwath,
                    ckd_spectral: tp.CKDSpectral) -> tp.CKDRadiometric:
    """Generate radiometric CKD."""
    rad_corr = np.empty((len(ckd_swath.act_angles),
                         len(ckd_spectral.wavelengths)))
    rad_corr[:] = 4222381848451.05
    return tp.CKDRadiometric(rad_corr)


def gen_pixel_mask(conf: dict,
                   ckd_dark: tp.CKDDark,
                   ckd_prnu: tp.CKDPRNU) -> npt.NDArray[np.bool_]:
    """Generate detector bad pixel mask."""
    mask = np.logical_or(abs(ckd_dark.current) > 1e4,
                         abs(ckd_dark.offset) > 600.0)
    mask = np.logical_or(mask, ckd_prnu.prnu_qe < 1e-3)
    return mask


def write_ckd(filename: str, ckd: tp.CKD, ckd_stray_file: str) -> None:
    """Write CKD contents to NetCDF file.

    Parameters
    ----------
    filename
        NetCDF file nam
    ckd
        CKD
    ckd_stray_file
        Stray light CKD is read and copied straight from a
        file. Converting the contents into a CKDStray object and back
        would be too much work.

    """
    nc = Dataset(filename, 'w')
    # Global attributes
    nc.title = "Tango Carbon calibration key data"
    nc.Conventions = "CF-1.11"
    nc.institution = "SRON Netherlands Institute for Space Research"
    nc.project = "Tango"
    nc.creator_name = "SRON/Earth Science"

    # Basic detector dimensions
    n_rows, n_cols = ckd.dark.offset.shape
    nc.createDimension('detector_row', n_rows)
    nc.createDimension('detector_column', n_cols)
    nc.createDimension('across_track_sample', len(ckd.swath.act_angles))
    nc.createDimension('wavelength', len(ckd.spectral.wavelengths))
    nc.createDimension('vector', 3)

    var = nc.createVariable('pixel_mask',
                            'b',
                            ('detector_row', 'detector_column'),
                            fill_value=np.int8(-128))
    var.long_name = 'detector pixel mask'
    var.valid_range = (np.int8(0), np.int8(1))
    var.flag_values = (np.int8(0), np.int8(1))
    var.flag_meanings = 'good bad'
    var[:] = ckd.pixel_mask

    # Dark
    grp = nc.createGroup('dark')
    var = grp.createVariable(
        'offset', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'detector offset'
    var.units = 'counts'
    var[:] = ckd.dark.offset
    var = grp.createVariable(
        'current', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'detector dark current'
    var.units = 'counts s-1'
    var[:] = ckd.dark.current

    # Noise
    grp = nc.createGroup('noise')
    var = grp.createVariable('g', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'noise conversion gain'
    var.description = 'noise model: sigma = sqrt(g*signal + n^2)'
    var.units = 'counts electrons-1'
    var[:] = ckd.noise.conversion_gain
    var = grp.createVariable('n', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'read noise'
    var.description = 'noise model: sigma = sqrt(g*signal + n^2)'
    var.units = 'counts^2 electrons-1'
    var[:] = ckd.noise.read_noise

    # Nonlinearity
    grp = nc.createGroup('nonlinearity')
    grp.createDimension('knots', len(ckd.nonlin.observed))
    var = grp.createVariable('observed', 'f8', 'knots')
    var.long_name = 'observed signal (spline knots)'
    var.units = 'counts'
    var[:] = ckd.nonlin.observed
    var = grp.createVariable('expected', 'f8', 'knots')
    var.long_name = 'expected (linear) signal (spline values)'
    var.units = 'counts'
    var[:] = ckd.nonlin.expected

    # PRNU
    grp = nc.createGroup('prnu')
    var = grp.createVariable('prnu', 'f8', ('detector_row', 'detector_column'))
    var.long_name = (
        'product of photoresponse non-uniformity and quantum efficiency')
    var.units = '1'
    var[:] = ckd.prnu.prnu_qe

    # Stray
    grp = nc.createGroup('stray')
    nc_in = Dataset(ckd_stray_file)
    grp.createDimension('kernel', nc_in['stray'].dimensions['kernel'].size)
    grp.createDimension('kernels_fft_size',
                        nc_in['stray'].dimensions['kernels_fft_size'].size)
    grp.createDimension('edges_of_box', 4)
    var = grp.createVariable('kernel_rows', 'i4', 'kernel')
    var.long_name = 'number of rows in each kernel'
    var[:] = nc_in['stray']['kernel_rows'][:].data
    var = grp.createVariable('kernel_cols', 'i4', 'kernel')
    var.long_name = 'number of cols in each kernel'
    var[:] = nc_in['stray']['kernel_cols'][:].data
    var = grp.createVariable('kernel_fft_sizes', 'i4', 'kernel')
    var.long_name = 'sizes of kernel FFTs'
    var.units = '16 bytes'
    var[:] = nc_in['stray']['kernel_fft_sizes'][:].data
    var = grp.createVariable('kernels_fft', 'f8', 'kernels_fft_size')
    var.long_name = 'Fourier transforms of the kernels'
    var[:] = nc_in['stray']['kernels_fft'][:].data
    var = grp.createVariable(
        'eta', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'internal scattering factor'
    var[:] = nc_in['stray']['eta'][:].data
    var = grp.createVariable(
        'weights', 'f8', ('kernel', 'detector_row', 'detector_column'))
    var.long_name = 'kernel weights'
    var[:] = nc_in['stray']['weights'][:].data
    var = grp.createVariable(
        'edges', 'i4', ('kernel', 'edges_of_box'))
    var.long_name = 'distances of subimage edges from the detector edges'
    var[:] = nc_in['stray']['edges'][:].data

    # Swath
    grp = nc.createGroup('swath')
    var = grp.createVariable('act_angle', 'f8', 'across_track_sample')
    var.long_name = 'across track angles'
    var.units = 'deg'
    var[:] = ckd.swath.act_angles
    var = grp.createVariable(
        'act_map', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'ACT angle of each detector pixel'
    var.units = 'deg'
    var[:] = ckd.swath.act_map
    var = grp.createVariable(
        'wavelength_map', 'f8', ('detector_row', 'detector_column'))
    var.long_name = 'wavelength of each detector pixel'
    var.units = 'nm'
    var[:] = ckd.swath.wavelength_map
    var = grp.createVariable(
        'row_map', 'f8', ('across_track_sample', 'wavelength'))
    var.long_name = 'row index of each L1B spectral element'
    var.units = '1'
    var[:] = ckd.swath.row_map
    var = grp.createVariable(
        'col_map', 'f8', ('across_track_sample', 'wavelength'))
    var.long_name = 'column index of each L1B spectral element'
    var.units = '1'
    var[:] = ckd.swath.col_map
    var = grp.createVariable(
        'line_of_sight', 'f8', ('across_track_sample', 'vector'))
    var.long_name = 'line of sight vector of each L1B spectrum'
    var.units = '1'
    var[:] = ckd.swath.line_of_sights

    # Spectral
    grp = nc.createGroup('spectral')
    var = grp.createVariable('wavelength', 'f8', 'wavelength')
    var.long_name = 'wavelengths of L1B spectra'
    var.units = 'nm'
    var[:] = ckd.spectral.wavelengths

    # Radiometric
    grp = nc.createGroup('radiometric')
    var = grp.createVariable(
        'radiometric', 'f8', ('across_track_sample', 'wavelength'))
    var.long_name = 'radiometric calibration constant'
    var[:] = ckd.radiometric.rad_corr

    nc.close()


def gen_ckd(config_user: dict | None = None) -> None:
    """Generate Tango Carbon CKD.

    Parameters
    ----------
    config_user
        Configuration dictionary directly from file, as given by the
        user, to be expanded with default values for parameters not
        specified by the user.

    """
    print_heading('Tango Carbon CKD generator', empty_line=False)
    print_system_info()
    print(flush=True)

    conf = merge_config_with_default(config_user, 'teds.ckd.ckd_generation')
    check_config(conf)

    print_heading('Generating CKD')

    log.info('Generating dark offset and current')
    ckd_dark = gen_dark(conf)

    log.info('Generating noise CKD')
    ckd_noise = gen_noise(conf)

    log.info('Generating nonlinearity CKD')
    ckd_nonlin = gen_nonlin(conf)

    log.info('Generating spectral CKD')
    ckd_spectral = gen_spectral(conf)

    log.info('Generating swath CKD')
    ckd_swath = gen_swath(conf, ckd_spectral)

    log.info('Generating PRNU and quantum efficiency CKD')
    ckd_prnu = gen_prnu(conf, ckd_swath, ckd_spectral)

    log.info('Generating radiometric CKD')
    ckd_radiometric = gen_radiometric(conf, ckd_swath, ckd_spectral)

    log.info('Generating detector bad pixel mask')
    mask = gen_pixel_mask(conf, ckd_dark, ckd_prnu)

    ckd_stray_dummy = tp.CKDStray([np.zeros((), dtype=np.complex128)],
                                  np.zeros(()),
                                  np.zeros(()),
                                  np.zeros((), dtype=np.int32))

    ckd = tp.CKD(ckd_dark.offset.shape[0],
                 ckd_dark.offset.shape[1],
                 mask,
                 ckd_dark,
                 ckd_noise,
                 ckd_nonlin,
                 ckd_prnu,
                 ckd_stray_dummy,
                 ckd_swath,
                 ckd_spectral,
                 ckd_radiometric)

    log.info('Writing output')
    write_ckd(conf['io_files']['ckd'], ckd, conf['io_files']['stray'])

    print_heading('Success')
