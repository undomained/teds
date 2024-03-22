"""
Wavelength Map for Spatial Samples
Based on Dispersion data
"""

import numpy as np
from tango.ckd_generation.generator_class import *
from netCDF4 import Dataset
from pandas import read_excel
import os
from scipy.interpolate import CubicSpline, interp1d

def import_data(cfg):
    filename = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_InstrumentModelInputs_20240307.xlsx'
    df = read_excel(filename, 'Dispersion', skiprows = 7, nrows=8,  usecols= 'B:G')
    # Tidying up data from excel file
    df.rename(columns = {"ACT field @ slit (mm)": "act_pos_at_slit"}, inplace = True)
    wl_cols = df.columns[1:]
    
    # Save data to netcdf file
    dispersion_nc = Dataset(cfg['paths']['dir_input']+'dispersion_data.nc', 'w', format="NETCDF4")
    dispersion_nc.createDimension('wavelengths', len(wl_cols))
    dispersion_nc.createDimension('act_pos_at_slit', len(df['act_pos_at_slit']))
    var_actpos = dispersion_nc.createVariable('act_pos_at_slit', 'float32', ['act_pos_at_slit'])
    var_actpos[:] = df['act_pos_at_slit']
    var_actpos.setncattr('units', 'mm')
    var_wavelength = dispersion_nc.createVariable('wavelengths', 'float32', ['wavelengths'])
    var_wavelength[:] = wl_cols
    var_actpos.setncattr('units', 'nm')
    var = dispersion_nc.createVariable('dispersion', 'float32', ['act_pos_at_slit','wavelengths'])
    var[:] = np.array(df[wl_cols])
    var.setncattr('units', "um/nm")
    var.setncattr('long_name', "nominal_linear_dispersion")
    dispersion_nc.close()


def calculate(cfg):
    filepath = cfg['paths']['dir_input']+'dispersion_data.nc'
    if not os.path.isfile(filepath):
        import_data(cfg)
    dispersion_nc = Dataset(filepath, 'a', format="NETCDF4")
    print('[wave_target] >> Calculating Dispersion Matrix')

    w = np.array(dispersion_nc['wavelengths'])
    y = np.array(dispersion_nc['act_pos_at_slit'])
    x = np.array(dispersion_nc['dispersion']) # position of projection
    # Take only positive values of act_pos because of symmetry
    x = x[y >= 0]
    y = y[y >= 0] 
    y = y / np.max(y)  # Normalize act_pos_at_slit to act_pos 
    # Normalize dispersion to image fraction
    x = x - np.min(x)
    x = x / np.max(x)
    
    # Get dimension of target 
    N_spat_px = cfg['dimensions']['spatial_samples_per_image']
    N_spec_px = cfg['dimensions']['spectral_detector_pixels']
    x_intp = np.linspace(0, 1, N_spec_px) # image fraction in spectral direction
    y_intp = np.abs(np.linspace(-1, 1, N_spat_px)) # image fraction in spatial direction, mirrored at nadir

    # Interpolate dispersion to number of spectral detector pixels, linear interpolation is used because the lines are fairly straight
    w_intp_x = np.zeros((len(y), N_spec_px)) # wavelength values to assign to the samples
    for i, xrow in enumerate(x):
        spl = interp1d(xrow, w, fill_value='extrapolate')
        w_intp_x[i] = spl(x_intp) # wavelength with act_pos (rows) and spec_pixels (cols)

    # Now interpolate ACT pos to number of samples with a Cubic Spline
    dispersion = np.zeros((N_spat_px, N_spec_px))
    for i, _ in enumerate(x_intp):
        spl = CubicSpline(y, w_intp_x[:,i], extrapolate=True)
        dispersion[:,i] = spl(y_intp)

    # Add to netcdf data file
    dispersion_nc.createGroup('interpolated')
    dispersion_nc['interpolated'].createDimension('spatial_samples_per_image', N_spat_px)
    dispersion_nc['interpolated'].createDimension('spectral_detector_pixels', N_spec_px)
    dispersion_var = dispersion_nc['interpolated'].createVariable('dispersion_interpolated', 'float32', ['spatial_samples_per_image', 'spectral_detector_pixels'])
    dispersion_var[:] = dispersion
    dispersion_var.setncattr('comment', 'minimum dispersion mapped to first col, max dispersion mapped to last col, must be changed later')




def generate(cfg):
    gen = ckd_generator()

    filepath = cfg['paths']['dir_input'] + 'dispersion_data.nc'
    if not os.path.isfile(filepath):
        calculate(cfg)
    wave_map_nc = Dataset(filepath, 'r', format="NETCDF4")
    try:
        wave_map = wave_map_nc['interpolated/dispersion_interpolated']
    except:
        # If data not found, recalculate
        wave_map_nc.close()
        calculate(cfg)
        wave_map_nc = Dataset(filepath, 'r', format="NETCDF4")
        wave_map = wave_map_nc['interpolated/dispersion_interpolated']

    dim_spat = 'spatial_samples_per_image'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    gen.dtype = 'float32'
    gen.data = wave_map
    gen.attr_names = ['comment', 'units']
    gen.attr_vals.append('minimum dispersion mapped to first col, max dispersion mapped to last col, must be changed later')
    gen.attr_vals.append('')
    return gen

