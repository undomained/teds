#!/bin/env python
# This source code is licensed under the 3-clause BSD license found
# in the LICENSE file in the root directory of this project.
'''CKD generator for Tango Carbon

This script compiles a calibration key data (CKD) file using
experimental or simulated calibration data as input. Depending on the
CKD step, the data could be copied directly from an input file (raw
CKD) or undergo some modification first such as interpolation to a
different temperature or unit conversion.

Usage: python ckd_generator.py ckd.yaml

'''

from netCDF4 import Dataset
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator
from teds import log
import numpy as np
import yaml
import sys

default_fill_value = -32767


def gen_header(conf: dict, nc_ckd: Dataset) -> None:
    """Generate basic CKD dimensions and global attributes.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    # Basic detector dimensions
    nc_ckd.createDimension('detector_row', conf['main']['n_detector_rows'])
    nc_ckd.createDimension('detector_column', conf['main']['n_detector_cols'])
    nc_ckd.createDimension('across_track_sample', conf['main']['n_act'])
    nc_ckd.createDimension('vector', 3)

    # Global attributes
    nc_ckd.title = "Tango Carbon calibration key data"
    nc_ckd.Conventions = "CF-1.11"
    nc_ckd.institution = "SRON Netherlands Institute for Space Research"
    nc_ckd.project = "Tango"
    nc_ckd.creator_name = "SRON/Earth Science"


def gen_dark(conf: dict, nc_ckd: Dataset) -> None:
    """Generate detector dark CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_grp = nc_ckd.createGroup('dark')
    nc_ckd_in = Dataset(conf['dark']['ckd_in'])
    nc_var = nc_grp.createVariable(
        'offset', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = 'detector offset'
    nc_var.units = 'counts'
    spline = PchipInterpolator(nc_ckd_in['temperature'][:],
                               np.nan_to_num(nc_ckd_in['offset'][:]),
                               0)
    nc_var[:] = spline(conf['main']['temperature'])
    nc_ckd_in = Dataset(conf['dark']['ckd_in'])
    nc_var = nc_grp.createVariable(
        'current', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = 'detector dark current'
    nc_var.units = 'counts s-1'
    spline = PchipInterpolator(nc_ckd_in['temperature'][:],
                               np.nan_to_num(nc_ckd_in['dark_current'][:]),
                               0)
    nc_var[:] = spline(conf['main']['temperature'])


def gen_noise(conf: dict, nc_ckd: Dataset) -> None:
    """Generate detector noise CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_grp = nc_ckd.createGroup('noise')
    nc_ckd_in = Dataset(conf['dark']['ckd_in'])
    nc_var = nc_grp.createVariable(
        'g', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = 'noise conversion gain'
    nc_var.description = 'noise model: sigma = sqrt(g*signal + n^2)'
    nc_var.units = 'counts electrons-1'
    nc_var[:] = (CubicSpline(nc_ckd_in['temperature'][:],
                             1 / nc_ckd_in['conversion_gain'][:])
                 (conf['main']['temperature']))
    nc_var = nc_grp.createVariable(
        'n', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = 'read noise'
    nc_var.description = 'noise model: sigma = sqrt(g*signal + n^2)'
    nc_var.units = 'counts^2 electrons-1'
    spline = PchipInterpolator(nc_ckd_in['temperature'][:],
                               np.nan_to_num(nc_ckd_in['read_noise'][:]),
                               0)
    nc_var[:] = spline(conf['main']['temperature'])


def gen_nonlin(conf: dict, nc_ckd: Dataset) -> None:
    """Generate detector nonlinearity CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_grp = nc_ckd.createGroup('nonlinearity')
    n_knots = 100
    nc_grp.createDimension('knots', n_knots)
    nc_var = nc_grp.createVariable('knots', 'f8', 'knots')
    nc_var.long_name = 'spline knots of the measured signal'
    nc_var.units = 'counts'
    nc_var[:] = np.linspace(-100, 10000, n_knots)
    nc_var = nc_grp.createVariable('y', 'f8', 'knots')
    nc_var.long_name = 'spline values of the measured signal'
    nc_var.units = 'counts'
    nc_var[:] = np.linspace(-100, 10000, n_knots)


def gen_swath_spectral(conf: dict, nc_ckd: Dataset) -> None:
    """Generate swath and spectral CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    act = (-1.7146, -1.1409, -0.56762, 0.00584267, 0.57931, 1.1526, 1.7262)
    wavelengths = (1.59, 1.61, 1.6325, 1.65, 1.675)
    wavelengths_new = (1.585, 1.61, 1.6325, 1.65, 1.679)
    row_indices_in = [
        [-3.760067, -3.755833, -3.749533, -3.746533, -3.7609],
        [-2.511, -2.508033, -2.5034, -2.500633, -2.508167],
        [-1.257133, -1.2556, -1.253167, -1.251533, -1.2547],
        [-7.68233e-05, -5.8018e-05, -3.853933e-05, -1.942e-05, 2.08867e-06],
        [1.257, 1.2555, 1.2531, 1.2515, 1.2547],
        [2.510867, 2.5079, 2.5033, 2.5006, 2.5082],
        [3.7599, 3.755733, 3.749467, 3.746533, 3.7609],
       ]
    col_indices_in = [
        [3.9143, 1.9941, 0.0393533, -1.961433, -4.54967],
        [4.041167, 2.1146, 0.1542133, -1.85, -4.441167],
        [4.1177, 2.1873, 0.22349, -1.7841, -4.3759],
        [4.1433, 2.2116, 0.2466633, -1.7618, -4.354133],
        [4.1178, 2.1873, 0.22352, -1.7841, -4.3759],
        [4.041267, 2.1147, 0.1542767, -1.850733, -4.4411],
        [3.9144, 1.9942, 0.03944433, -1.961367, -4.549567],
       ]
    pixel_size = 0.015  # mm
    row_indices = np.asarray(row_indices_in) / pixel_size
    col_indices = np.asarray(col_indices_in) / pixel_size
    # The spot measurement data assumes a 512x640 detector. For any
    # other dimensions the generated pixel indices need to be scaled
    # appropriately.
    row_indices *= conf['main']['n_detector_rows'] / 512
    col_indices *= conf['main']['n_detector_cols'] / 640
    # Shift to the detector center
    row_indices += conf['main']['n_detector_rows'] / 2
    col_indices += conf['main']['n_detector_cols'] / 2

    for i_act in range(len(act)):
        row_indices[i_act, :] = (
            CubicSpline(wavelengths, row_indices[i_act, :])(wavelengths_new))
        col_indices[i_act, :] = (
            CubicSpline(wavelengths, col_indices[i_act, :])(wavelengths_new))
    wavelengths_new = wavelengths
    los_spline = bisplrep(np.asarray(5 * [act]).transpose().flatten(),
                          col_indices.flatten(),
                          row_indices.flatten(),
                          kx=3,
                          ky=3)
    min_act = -1.714
    max_act = 1.714
    act_angles = np.linspace(
        min_act, max_act, nc_ckd.dimensions['across_track_sample'].size)
    row_indices = bisplev(act_angles, np.arange(0.0, 640.0, 1.0), los_spline)
    nc_grp = nc_ckd.createGroup('swath')
    nc_var = nc_grp.createVariable(
        'row_index', 'f8', ('across_track_sample', 'detector_column'))
    nc_var.long_name = 'row indices at which to extract the spectra'
    nc_var[:] = np.flip(np.flip(row_indices, 0), 1)
    nc_var = nc_grp.createVariable(
        'line_of_sight', 'f8', ('across_track_sample', 'vector'))
    nc_var[:] = (1.0, 0.0, 0.0)
    nc_var.long_name = 'line of sight vectors'
    # Spectral
    wave_spline = bisplrep(np.asarray(5 * [act]).transpose().flatten(),
                           col_indices.flatten(),
                           np.asarray(7 * [wavelengths]).flatten(),
                           kx=3,
                           ky=3)
    wavelength = bisplev(act_angles, np.arange(0.0, 640.0, 1.0), wave_spline)
    nc_grp = nc_ckd.createGroup('spectral')
    nc_var = nc_grp.createVariable(
        'wavelength', 'f8', ('across_track_sample', 'detector_column'))
    nc_var.long_name = 'wavelengths of L1B spectra'
    nc_var.units = 'nm'
    wavelength *= 1e3
    nc_var[:] = wavelength


def gen_prnu(conf: dict, nc_ckd: Dataset) -> None:
    """Generate PRNU x QE CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_ckd_in = Dataset(conf['prnu']['ckd_in'])
    # Number of wavelengths for quantum efficiency (QE)
    n_waves = nc_ckd_in.dimensions['wavelength'].size
    # Number of L1B spectra
    n_act = nc_ckd.dimensions['across_track_sample'].size
    # QE is more difficult because it is provided per wavelength.
    # First step is to interpolate QE onto the target temperature.
    qe = np.zeros(n_waves)
    temperatures_qe = nc_ckd_in['temperature_qe'][:]
    for i in range(n_waves):
        spline = CubicSpline(temperatures_qe,
                             nc_ckd_in['quantum_efficiency'][:, i])
        qe[i] = spline(conf['main']['temperature'])
    # Next interpolate QE onto the L1B spectra
    qe_spline = CubicSpline(nc_ckd_in['wavelength'][:], qe)
    spectra = np.zeros((n_act, conf['main']['n_detector_cols']))
    wavelength = nc_ckd['spectral/wavelength'][:]
    for i_act in range(n_act):
        spectra[i_act, :] = qe_spline(wavelength[i_act, :])
    # At this point, we have a QE value for each spectrum (FOV) index and
    # each detector column. Final step is to interpolate, per column, the
    # QE values from FOV indices to detector rows using the FOV CKD.
    # QE values of all pixels
    qe_map = np.zeros((conf['main']['n_detector_rows'],
                       conf['main']['n_detector_cols']))
    # Interpolation target is a list of integer row numbers
    row_indices_out = np.arange(0.0, conf['main']['n_detector_rows'], 1.0)
    row_indices = nc_ckd['swath/row_index'][:]
    row_indices = np.flip(row_indices, 0)
    for i_col in range(conf['main']['n_detector_cols']):
        drawing_spline = CubicSpline(row_indices[:, i_col],
                                     spectra[:, i_col])
        qe_map[:, i_col] = drawing_spline(row_indices_out)
    spline = PchipInterpolator(nc_ckd_in['temperature'][:],
                               np.nan_to_num(nc_ckd_in['prnu'][:]),
                               0)
    prnu = spline(conf['main']['temperature'])
    # Combine PRNU and QE into one variable (for now)
    prnu = prnu * qe_map
    nc_grp = nc_ckd.createGroup('prnu')
    nc_var = nc_grp.createVariable(
        'prnu', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = (
        'product of photoresponse non-uniformity and quantum efficiency')
    nc_var.units = '1'
    nc_var[:] = prnu


def gen_stray(conf: dict, nc_ckd: Dataset) -> None:
    """Generate stray light CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_grp = nc_ckd.createGroup('stray')
    nc_ckd_in = Dataset(conf['stray']['ckd_in'])['stray']

    nc_grp.createDimension('kernel', nc_ckd_in.dimensions['kernel'].size)
    nc_grp.createDimension('kernels_fft_size',
                           nc_ckd_in.dimensions['kernels_fft_size'].size)
    nc_grp.createDimension('edges_of_box', 4)

    nc_var = nc_grp.createVariable('kernel_rows', 'i4', 'kernel')
    nc_var.long_name = 'number of rows in each kernel'
    nc_var[:] = nc_ckd_in['kernel_rows'][:]

    nc_var = nc_grp.createVariable('kernel_cols', 'i4', 'kernel')
    nc_var.long_name = 'number of cols in each kernel'
    nc_var[:] = nc_ckd_in['kernel_cols'][:]

    nc_var = nc_grp.createVariable('kernel_fft_sizes', 'i4', 'kernel')
    nc_var.long_name = 'sizes of kernel FFTs'
    nc_var.units = '16 bytes'
    nc_var[:] = nc_ckd_in['kernel_fft_sizes'][:]

    nc_var = nc_grp.createVariable('kernels_fft', 'f8', 'kernels_fft_size')
    nc_var.long_name = 'Fourier transforms of the kernels'
    nc_var[:] = nc_ckd_in['kernels_fft'][:]

    nc_var = nc_grp.createVariable(
        'eta', 'f8', ('detector_row', 'detector_column'))
    nc_var.long_name = 'internal scattering factor'
    nc_var[:] = nc_ckd_in['eta'][:]

    nc_var = nc_grp.createVariable(
        'weights', 'f8', ('kernel', 'detector_row', 'detector_column'))
    nc_var.long_name = 'kernel weights'
    nc_var[:] = nc_ckd_in['weights'][:]

    nc_var = nc_grp.createVariable(
        'edges', 'f8', ('kernel', 'edges_of_box'))
    nc_var.long_name = 'distances of subimage edges from the detector edges'
    nc_var[:] = nc_ckd_in['edges'][:]


def gen_radiometric(conf: dict, nc_ckd: Dataset) -> None:
    """Generate radiometric CKD and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_grp = nc_ckd.createGroup('radiometric')
    nc_var = nc_grp.createVariable(
        'radiometric', 'f8', ('across_track_sample', 'detector_column'))
    nc_var.long_name = 'radiometric calibration constant'
    nc_var[:] = 4222381848451.05


def gen_pixel_mask(conf: dict, nc_ckd: Dataset) -> None:
    """Generate detector bad pixel mask and save into a NetCDF file.

    Parameters
    ----------
    conf
        settings read from the YAML file
    nc_ckd
        NetCDF file to save the result to

    """
    nc_var = nc_ckd.createVariable('pixel_mask',
                                   'b',
                                   ('detector_row', 'detector_column'),
                                   fill_value=np.int8(-128))
    nc_var.long_name = 'detector pixel mask'
    nc_var.valid_range = (np.int8(0), np.int8(1))
    nc_var.flag_values = (np.int8(0), np.int8(1))
    nc_var.flag_meanings = 'good bad'
    mask = np.logical_or(abs(nc_ckd['dark/current'][:]) > 1e4,
                         abs(nc_ckd['dark/offset'][:]) > 600.0)
    mask = np.logical_or(mask, abs(nc_ckd['prnu/prnu'][:]) < 1e-10)
    nc_var[:] = mask


def gen_ckd(conf: dict) -> None:
    print('##############################\n'
          '# Tango Carbon CKD generator #\n'
          '##############################')

    nc_ckd = Dataset(conf['main']['ckd_out'], 'w')

    log.info('Generating dimensions and global attributes')
    gen_header(conf, nc_ckd)

    print('Generating dark offset and current')
    gen_dark(conf, nc_ckd)

    print('Generating noise CKD')
    gen_noise(conf, nc_ckd)

    print('Generating nonlinearity CKD')
    gen_nonlin(conf, nc_ckd)

    print('Generating swath and spectral CKD')
    gen_swath_spectral(conf, nc_ckd)

    print('Generating PRNU and quantum efficiency CKD')
    gen_prnu(conf, nc_ckd)

    print('Generating stray light CKD')
    gen_stray(conf, nc_ckd)

    print('Generating radiometric CKD')
    gen_radiometric(conf, nc_ckd)

    print('Generating detector bad pixel mask')
    gen_pixel_mask(conf, nc_ckd)

    nc_ckd.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage: ckd_generator.py [ckd.yaml]')
        exit(0)

    conf = yaml.safe_load(open(sys.argv[1]))
    gen_ckd(conf)
