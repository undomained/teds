#!/bin/env python

from netCDF4 import Dataset
from scipy.interpolate import CubicSpline
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.optimize import curve_fit
import numpy as np
import yaml
import os
import sys

if len(sys.argv) < 2:
    print('usage: ckd_generator.py [ckd.yaml]')
    exit(0)

conf = yaml.safe_load(open(sys.argv[1]))

if conf['main']['append_to_ckd']:
    nc_ckd = Dataset(conf['main']['ckd_out'], 'r+')
else:
    nc_ckd = Dataset(conf['main']['ckd_out'], 'w')

# Basic detector dimensions
n_rows = conf['main']['n_detector_rows']
n_cols = conf['main']['n_detector_cols']
detector_shape = (n_rows, n_cols)
nc_ckd.createDimension('spatial_detector_pixels', n_rows)
nc_ckd.createDimension('spectral_detector_pixels', n_cols)
nc_ckd.createDimension('vector_elements', 3)
nc_ckd.createDimension('number_of_views', 1)
nc_ckd.createDimension('spatial_samples_per_image', 100)

# Global attributes
nc_ckd.Conventions = "CF-1.6"
nc_ckd.institution = "SRON Netherlands Institute for Space Research"
nc_ckd.instrument = "Tango"
nc_ckd.project = "Tango"
nc_ckd.creator_name = "SRON/Earth Science"


def interp(temperatures, data_in):
    '''Interpolate detector CKD for target temperature

    Arguments:
      data_in - a temperature dependent detector CKD (n_temps, n_rows, n_cols)

    Returns:
      CKD (n_rows, n_cols) interpolated at the target temperature

    '''
    target_temperature = conf['main']['temperature']
    data_out = np.zeros(detector_shape)
    for i in range(n_rows):
        for j in range(n_cols):
            if data_in[:, i, j].mask.any():
                continue
            data_out[i, j] = CubicSpline(
                temperatures, data_in[:, i, j])(target_temperature)
    return data_out


def interpExp(temperatures, data_in):
    '''Fit exponential and interpolate to target temperature

    If the delta_T array does not exist, compute and save it to disk
    first. Otherwise read it from disk.

    Arguments:
      data_in - a detector CKD with the dimensions (n_temps, n_rows, n_cols)

    Returns:
      CKD (n_rows, n_cols) interpolated at the target temperature

    '''
    delta_T_filename = 'delta_T.bin'
    if os.path.isfile(delta_T_filename):
        delta_T = np.fromfile(delta_T_filename)
        delta_T = delta_T.reshape(detector_shape)
    else:
        delta_T = np.zeros(detector_shape)
        for i in range(n_rows):
            print('Fitting delta_T for detector row', i)
            for j in range(n_cols):
                # Omit the first temperature
                x_data = temperatures[1:]
                y_data = data_in[1:, i, j]
                # Skip masked pixels
                if y_data.mask.any():
                    continue
                I_ref = data_in[-1, i, j]
                T_ref = 15.0

                def exponential(T, delta_T):
                    return I_ref * 2**((T - T_ref) / delta_T)

                pars, cov = curve_fit(exponential, x_data, y_data)
                delta_T[i, j] = pars[0]
        delta_T.tofile('delta_T.bin')

    data_out = np.zeros(detector_shape)
    for i in range(n_rows):
        for j in range(n_cols):
            # Skip bad pixels
            if data_in[:, i, j].mask.any():
                continue
            I_ref = data_in[-1, i, j]

            def exponential(T, I_ref, T_ref, delta_T):
                return I_ref * 2**((T - T_ref) / delta_T)

            data_out[i, j] = exponential(conf['temperature'],
                                         data_in[-1, i, j],
                                         15.0,
                                         delta_T[i, j])
    return data_out


# Pixel mask
mask = np.zeros(detector_shape)
nc_mask = nc_ckd.createVariable(
    'mask', 'u1', ('spatial_detector_pixels', 'spectral_detector_pixels'))
nc_mask.long_name = 'detector pixel mask'
nc_mask.flag_values = (np.uint8(0), np.uint8(1))
nc_mask.flag_meanings = 'good bad'

nc_vp_mask = nc_ckd.createVariable(
    'vp_mask', 'u1', 'number_of_views')
nc_vp_mask.comment = 'unused'

# Create all groups. These are optionally filled later.
nc_ckd.createGroup('DARK')
nc_ckd.createGroup('NON_LINEARITY')
nc_ckd.createGroup('NOISE')
nc_ckd.createGroup('NON_LINEARITY')
nc_ckd.createGroup('PRNU')
nc_ckd.createGroup('STRAYLIGHT')
nc_ckd.createGroup('FIELD_OF_VIEW')
nc_ckd.createGroup('SWATH')
nc_ckd.createGroup('WAVELENGTH')
nc_ckd.createGroup('RADIOMETRIC')

nc_ckd['DARK'].createDimension('dark_number_of_coefficients', 1)

nc_ckd['DARK'].createVariable('dark_skip', 'u1')[:] = 0
nc_ckd['NOISE'].createVariable('noise_skip', 'u1')[:] = 0
nc_ckd['NON_LINEARITY'].createVariable('nonlin_skip', 'u1')[:] = 1
nc_ckd['PRNU'].createVariable('prnu_skip', 'u1')[:] = 0
nc_ckd['STRAYLIGHT'].createVariable('stray_skip', 'u1')[:] = 1
nc_ckd['SWATH'].createVariable('swath_skip', 'u1')[:] = 1
nc_ckd['RADIOMETRIC'].createVariable('rad_skip', 'u1')[:] = 0

# Dark offset
if 'dark_offset' in conf:
    print('Processing dark offset')
    nc_ckd_in = Dataset(conf['dark_offset']['ckd_in'])
    nc_offset = nc_ckd['DARK'].createVariable(
        'dark_offset',
        'f8',
        ('spatial_detector_pixels',
         'spectral_detector_pixels',
         'dark_number_of_coefficients'))
    nc_offset.long_name = 'detector offset'
    nc_offset.units = 'counts'
    nc_ckd['DARK/dark_offset'][:, :, 0] = interp(nc_ckd_in['temperature'][:],
                                                 nc_ckd_in['offset'][:])

# Dark current
if 'dark_current' in conf:
    print('Processing dark current')
    nc_ckd_in = Dataset(conf['dark_current']['ckd_in'])
    nc_offset = nc_ckd['DARK'].createVariable(
        'dark_current',
        'f8',
        ('spatial_detector_pixels',
         'spectral_detector_pixels',
         'dark_number_of_coefficients'))
    nc_offset.long_name = 'detector offset'
    nc_offset.units = 'counts'
    nc_ckd['DARK/dark_current'][:, :, 0] = interp(nc_ckd_in['temperature'][:],
                                                  nc_ckd_in['dark_current'][:])

# Dark current temperature
nc_ckd['DARK'].createDimension('double_single', 1)
nc_dark_nominal_temperature = nc_ckd['DARK'].createVariable(
    'dark_nominal_temperature',
    'f8',
    ('double_single'))
nc_dark_nominal_temperature.long_name = 'temperature'
nc_dark_nominal_temperature.units = 'degC'
nc_dark_nominal_temperature[0] = conf['main']['temperature']

# Noise
if 'dark_offset' in conf:
    print('Processing noise')
    nc_noise_g = nc_ckd['NOISE'].createVariable(
        'noise_g',
        'f8',
        ('spatial_detector_pixels', 'spectral_detector_pixels'))
    nc_noise_g.long_name = 'noise conversion gain'
    nc_noise_g.comment = 'noise model: sigma = sqrt(g*signal + n)'
    nc_noise_g.units = 'counts/e'
    nc_noise_g[:] = (CubicSpline(nc_ckd_in['temperature'][:],
                                 1 / nc_ckd_in['conversion_gain'][:])
                     (conf['main']['temperature']))

    nc_noise_n = nc_ckd['NOISE'].createVariable(
        'noise_n',
        'f8',
        ('spatial_detector_pixels', 'spectral_detector_pixels'))
    nc_noise_n.long_name = 'scaled read noise'
    nc_noise_n.comment = 'noise model: sigma = sqrt(g*signal + n)'
    nc_noise_n.units = 'counts^2/e'
    nc_noise_n[:] = interp(nc_ckd_in['temperature'][:],
                           nc_ckd_in['read_noise'][:])


# Field of view
print('Processing field of view CKD')
act = (-1.7146, -1.1409, -0.56762, 0.00584267, 0.57931, 1.1526, 1.7262)
wavelengths = (1.59, 1.61, 1.6325, 1.65, 1.675)
wavelengths_new = (1.585, 1.61, 1.6325, 1.65, 1.679)
row_indices = [
    [-3.760067, -3.755833, -3.749533, -3.746533, -3.7609],
    [-2.511, -2.508033, -2.5034, -2.500633, -2.508167],
    [-1.257133, -1.2556, -1.253167, -1.251533, -1.2547],
    [-7.68233e-05, -5.8018e-05, -3.853933e-05, -1.942e-05, 2.08867e-06],
    [1.257, 1.2555, 1.2531, 1.2515, 1.2547],
    [2.510867, 2.5079, 2.5033, 2.5006, 2.5082],
    [3.7599, 3.755733, 3.749467, 3.746533, 3.7609],
]
col_indices = [
    [3.9143, 1.9941, 0.0393533, -1.961433, -4.54967],
    [4.041167, 2.1146, 0.1542133, -1.85, -4.441167],
    [4.1177, 2.1873, 0.22349, -1.7841, -4.3759],
    [4.1433, 2.2116, 0.2466633, -1.7618, -4.354133],
    [4.1178, 2.1873, 0.22352, -1.7841, -4.3759],
    [4.041267, 2.1147, 0.1542767, -1.850733, -4.4411],
    [3.9144, 1.9942, 0.03944433, -1.961367, -4.549567],
]
pixel_size = 0.015  # mm
row_indices = np.asarray(row_indices) / pixel_size
col_indices = np.asarray(col_indices) / pixel_size
row_indices += conf['main']['n_detector_rows'] / 2
col_indices += conf['main']['n_detector_cols'] / 2

for i_act in range(len(act)):
    row_indices[i_act, :] = (
        CubicSpline(wavelengths, row_indices[i_act, :])(wavelengths_new))
    col_indices[i_act, :] = (
        CubicSpline(wavelengths, col_indices[i_act, :])(wavelengths_new))
wavelengths_new = wavelengths

fov_spline = bisplrep(np.asarray(5 * [act]).transpose().flatten(),
                      col_indices.flatten(),
                      row_indices.flatten(),
                      kx=3,
                      ky=3)

n_act = 100
min_act = -1.714
max_act = 1.714
fov_act_angles = np.linspace(min_act, max_act, n_act)
fov_ispat = bisplev(fov_act_angles, np.arange(0.0, 640.0, 1.0), fov_spline)

nc_fov_act_angles = nc_ckd['FIELD_OF_VIEW'].createVariable(
        'fov_act_angles',
        'f8',
        ('spatial_samples_per_image',))
nc_fov_act_angles.long_name = (
    'across track rotation stage position per spatial sample')
nc_fov_act_angles.units = 'radians'
nc_fov_act_angles[:] = fov_act_angles

nc_fov_ispat = nc_ckd['FIELD_OF_VIEW'].createVariable(
        'fov_ispat',
        'f8',
        ('spatial_samples_per_image', 'spectral_detector_pixels'))
nc_fov_ispat.long_name = (
    'floating point spatial detector pixel indices of the spectra')
nc_fov_ispat.units = '1'
nc_fov_ispat[:] = np.flip(np.flip(fov_ispat, 0), 1)

# Spectral
print('Processing spectral CKD')
wave_spline = bisplrep(np.asarray(5 * [act]).transpose().flatten(),
                       col_indices.flatten(),
                       np.asarray(7 * [wavelengths]).flatten(),
                       kx=3,
                       ky=3)
wave_target = bisplev(fov_act_angles, np.arange(0.0, 640.0, 1.0), wave_spline)

nc_wave_target = nc_ckd['WAVELENGTH'].createVariable(
        'wave_target',
        'f8',
        ('spatial_samples_per_image', 'spectral_detector_pixels'))
nc_wave_target.long_name = 'wavelengths of L1B spectra'
nc_wave_target.units = 'nm'
nc_wave_target[:] = wave_target * 1e3

nc_wave_spectra = nc_ckd['WAVELENGTH'].createVariable(
        'wave_spectra',
        'f8',
        ('spatial_samples_per_image', 'spectral_detector_pixels'))
nc_wave_spectra.comment = 'unused'

# PRNU
print('Processing PRNU')
# Number of wavelengths for quantum efficiency (QE)
n_waves = nc_ckd_in.dimensions['wavelength'].size
# Number of L1B spectra
n_act = nc_ckd.dimensions['spatial_samples_per_image'].size
# Temperatures at which the input CKD was measured (QE below has its
# own temperature list)
temperatures = nc_ckd_in['temperature'][:]
# QE is more difficult because it is provided per wavelength.
# First step is to interpolate QE onto the target temperature.
print('Interpolating quantum efficiency')
qe = np.zeros(n_waves)
temperatures_qe = nc_ckd_in['temperature_qe'][:]
for i in range(n_waves):
    spline = CubicSpline(temperatures_qe,
                         nc_ckd_in['quantum_efficiency'][:, i])
    qe[i] = spline(conf['main']['temperature'])
# Next interpolate QE onto the L1B spectra
qe_spline = CubicSpline(nc_ckd_in['wavelength'][:], qe)
spectra = np.zeros((n_act, n_cols))
for i_fov in range(n_act):
    spectra[i_fov, :] = qe_spline(wave_target[i_fov, :])
# At this point, we have a QE value for each spectrum (FOV) index and
# each detector column. Final step is to interpolate, per column, the
# QE values from FOV indices to detector rows using the FOV CKD.
# QE values of all pixels
qe_map = np.zeros(detector_shape)
# Interpolation target is a list of integer row numbers
row_indices_out = np.arange(0.0, 512.0, 1.0)
for i_col in range(n_cols):
    drawing_spline = CubicSpline(fov_ispat[:, i_col],
                                 spectra[:, i_col])
    qe_map[:, i_col] = drawing_spline(row_indices_out)
prnu = interp(nc_ckd_in['temperature'][:], nc_ckd_in['prnu'][:])
# Combine PRNU and QE into one variable (for now)
prnu = prnu * qe_map

nc_prnu = nc_ckd['PRNU'].createVariable(
    'prnu_prnu',
    'f8',
    ('spatial_detector_pixels', 'spectral_detector_pixels'))
nc_prnu.long_name = (
    'product of pixel response non-uniformity and quantum efficiency')
nc_prnu.units = '1'
nc_prnu[:] = prnu

# Radiometric
print('Processing radiometric CKD')
nc_rad = nc_ckd['RADIOMETRIC'].createVariable(
    'rad_spectra',
    'f8',
    ('spatial_samples_per_image', 'spectral_detector_pixels'))
nc_rad.long_name = 'radiometric calibration constants'
nc_rad.units = '1'
nc_rad[:] = 4222381848451.05

nc_ckd.close()
