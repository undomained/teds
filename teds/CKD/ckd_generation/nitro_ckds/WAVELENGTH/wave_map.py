"""
Wavelength Map with Spectral Smile
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *
import pandas as pd
from netCDF4 import Dataset
import os
from scipy.interpolate import interp1d


def import_data(cfg):
    filename = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_InstrumentModelInputs_20240307.xlsx'
    df = pd.read_excel(filename, 'SpectralSmile', skiprows = 6, nrows=21,  usecols= 'B:M')
    
    # Trying to import horrible excel file...
    obsolete_columns = [col for col in df.columns if 'delta' in str(col)]
    obsolete_columns.append('Centroid')
    df.drop(obsolete_columns, axis = 'columns', inplace=True)
    wl_cols = [col for col in df.columns if 'ACT' not in str(col)]
    # Values in dataframe assumed to be distance on detector in mm

    # Save data to netcdf file
    wave_map_nc = Dataset(cfg['paths']['dir_input']+'spectral_smile_data.nc', 'w', format="NETCDF4")
    wave_map_nc.createDimension('wavelengths', len(wl_cols))
    wave_map_nc.createDimension('act_pos', len(df['ACT Pos']))
    var_actpos = wave_map_nc.createVariable('act_pos', 'float32', ['act_pos'])
    var_actpos[:] = df['ACT Pos']
    var_wavelength = wave_map_nc.createVariable('wavelengths', 'float32', ['wavelengths'])
    var_wavelength[:] = wl_cols
    var = wave_map_nc.createVariable('shift', 'float32', ['act_pos','wavelengths'])
    var[:] = np.array(df[wl_cols])
    var.setncattr('units', "mm?")
    wave_map_nc.close()


def calculate(cfg):
    filepath = cfg['paths']['dir_input']+'spectral_smile_data.nc'
    if not os.path.isfile(filepath):
        import_data(cfg)
    wave_map_nc = Dataset(filepath, 'a', format="NETCDF4")

    print('[wave_map] >> Calculating Spectral Smile')

    w = np.array(wave_map_nc['wavelengths'])
    y = np.array(wave_map_nc['act_pos'])
    x = np.array(wave_map_nc['shift'])
    x = x / (cfg['image_extent_spec']/2) # Normalize to image fraction

    N_spat_px = cfg['dimensions']['spatial_detector_pixels']
    N_spec_px = cfg['dimensions']['spectral_detector_pixels']
    
    x_intp = np.linspace(-1, 1, N_spec_px) # detector pixels in spectral direction
    y_intp = np.abs(np.linspace(-1, 1, N_spat_px)) # detector pixels in spatial direction
    # the data are horizontally mirrored with the assumption that the map is symmetrical...
    w_intp_x = np.zeros((len(y), N_spec_px)) # wavelength values to assign to the pixels

    # First interpolate the wavelengths, linear interpolation is used
    for i, xrow in enumerate(x):
        spl = interp1d(xrow, w, fill_value='extrapolate')
        w_intp_x[i] = spl(x_intp)

    # Now interpolate over ACT pos
    wave_map = np.zeros((N_spat_px, N_spec_px))
    for i, _ in enumerate(x_intp):
        spl = interp1d(y, w_intp_x[:,i], fill_value='extrapolate')
        wave_map[:,i] = spl(y_intp)

    # Add to netcdf data file
    wave_map_nc.createGroup('interpolated')
    wave_map_nc['interpolated'].createDimension('spatial_detector_pixels', N_spat_px)
    wave_map_nc['interpolated'].createDimension('spectral_detector_pixels', N_spec_px)
    wave_map_var = wave_map_nc['interpolated'].createVariable('wave_map', 'float32', ['spatial_detector_pixels', 'spectral_detector_pixels'])
    wave_map_var[:] = wave_map


def generate(cfg):
    gen = ckd_generator()
    filepath = cfg['paths']['dir_input'] + 'spectral_smile_data.nc'
    if not os.path.isfile(filepath):
        calculate(cfg)
    wave_map_nc = Dataset(filepath, 'r', format="NETCDF4")
    try:
        wave_map = wave_map_nc['interpolated/wave_map']
    except:
        # If data not found, recalculate
        wave_map_nc.close()
        calculate(cfg)
        wave_map_nc = Dataset(filepath, 'r', format="NETCDF4")
        wave_map = wave_map_nc['interpolated/wave_map']

    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    gen.dtype = 'float32'
    gen.data = wave_map
    gen.attr_names = ['comment', 'units']
    gen.attr_vals = ['origin chosen such that center row (nadir) at 447.5nm has no shift', 'nm']

    return gen

