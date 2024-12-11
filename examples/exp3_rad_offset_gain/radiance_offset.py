import numpy as np
import sys
import yaml
import netCDF4 as nc
from copy import deepcopy
from teds.lib.libWrite import writevariablefromname

#   main program ##############################################################

def get_l1b(l1b_filename):
    # getting l1b data from file
    nc_l1b = nc.Dataset(l1b_filename, mode='r')
        
    l1b_data = {}
    l1b_data['wavelength'] = deepcopy(nc_l1b['observation_data/wavelength'][:])
    l1b_data['radiance'] = deepcopy(nc_l1b['observation_data/radiance'][:])
    l1b_data['noise'] = deepcopy(nc_l1b['observation_data/radiance_stdev'][:])
    l1b_data['sza'] = deepcopy(nc_l1b['geolocation_data/solar_zenith'][:])
    l1b_data['saa'] = deepcopy(nc_l1b['geolocation_data/solar_azimuth'][:])
    l1b_data['vza'] = deepcopy(nc_l1b['geolocation_data/sensor_zenith'][:])
    l1b_data['vaa'] = deepcopy(nc_l1b['geolocation_data/sensor_azimuth'][:])
    l1b_data['latitude'] = deepcopy(nc_l1b['geolocation_data/latitude'][:])
    l1b_data['longitude'] = deepcopy(nc_l1b['geolocation_data/longitude'][:])
    # Extract mask
    nc_var = nc_l1b['observation_data/radiance']
    l1b_data['mask'] = ~nc_var[:].mask
    if not nc_var[:].mask.shape:
        l1b_data['mask'] = np.full(nc_var[:].shape, True)
    
    return (l1b_data)

def sim_modified_output(filename, l1b_output):
    output = nc.Dataset(filename, mode='w')
    output.title = 'Tango Carbon level 1B data'

    nalt, nact, nwav = l1b_output['radiance'].shape
    output.createDimension('wavelength', nwav)
    output.createDimension('across_track_sample', nact)
    output.createDimension('along_track_sample', nalt)

    nc_grp = output.createGroup('geolocation_data')
    _dim2 = ('along_track_sample', 'across_track_sample')

    nc_var = nc_grp.createVariable(
        'latitude', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'latitude at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees_north'
    nc_var[:] = l1b_output['latitude']

    nc_var = nc_grp.createVariable(
        'longitude', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'longitude at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees_east'
    nc_var[:] = l1b_output['latitude']

    nc_var = nc_grp.createVariable(
        'solar_zenith', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'solar zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_output['sza']

    nc_var = nc_grp.createVariable(
        'solar_azimuth', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'solar azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_output['saa']

    nc_var = nc_grp.createVariable(
        'sensor_zenith', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'sensor zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_output['vza']

    nc_var = nc_grp.createVariable(
        'sensor_azimuth', 'f4', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'sensor azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_output['vaa']

    nc_grp = output.createGroup('observation_data')
    _dim3 = ('along_track_sample', 'across_track_sample', 'wavelength')
    _dim2 = ('across_track_sample', 'wavelength')
    nc_var = nc_grp.createVariable(
        'wavelength', 'f8', _dim2, fill_value=-32767.0)
    nc_var.long_name = 'wavelength'
    nc_var.valid_min = 0
    nc_var.valid_max = 8.e3
    nc_var.units = 'nm'
    nc_var[:] = l1b_output['wavelength']

    nc_var = nc_grp.createVariable(
        'radiance', 'f8', _dim3, fill_value=-32767.0)
    nc_var.long_name = 'spectral radiance'
    nc_var.valid_min = 0
    nc_var.valid_max = 1.e28
    nc_var.units = 'photons / (sr nm m2 s)'
    nc_var[:] = l1b_output['radiance']

    nc_var = nc_grp.createVariable(
        'radiance_stdev', 'f8', _dim3, fill_value=-32767.0)
    nc_var.long_name = 'relative radiance noise (1-sigma)'
    nc_var.units = 'photons/(sr nm m2 s)'
    nc_var.valid_min = 0.0
    nc_var.valid_max = 1e24
    nc_var[:] = l1b_output['noise']

    output.close()

def add_radiance_offset(filename_in, filename_out, rad_offset):

    l1b = get_l1b(filename_in)
    nalt, nact, nwav = l1b['radiance'].shape

    l1b_scaled = {}
    l1b_scaled['sza'] = np.empty([nalt, nact])
    l1b_scaled['saa'] = np.empty([nalt, nact])
    l1b_scaled['vza'] = np.empty([nalt, nact])
    l1b_scaled['vaa'] = np.empty([nalt, nact])
    l1b_scaled['latitude'] = np.empty([nalt, nact])
    l1b_scaled['longitude'] = np.empty([nalt, nact])
    l1b_scaled['noise'] = np.empty([nalt, nact, nwav])
    l1b_scaled['radiance'] = np.empty([nalt, nact, nwav])
    l1b_scaled['wavelength'] = np.empty([nwav])

    for ialt in range(nalt):
        for iact in range(nact):
            offset = rad_offset*np.max(l1b['radiance'][0, iact, :])
            l1b_scaled['radiance'][ialt, iact, :] = l1b['radiance'][0, iact, :]+offset
        l1b_scaled['sza'][ialt, :] = l1b['sza'][0, :]
        l1b_scaled['saa'][ialt, :] = l1b['saa'][0, :]
        l1b_scaled['vza'][ialt, :] = l1b['vza'][0, :]
        l1b_scaled['vaa'][ialt, :] = l1b['vaa'][0, :]
        l1b_scaled['latitude'][ialt, :] = l1b['latitude'][0, :]
        l1b_scaled['longitude'][ialt, :] = l1b['longitude'][0, :]
        l1b_scaled['noise'][ialt, :, :] = l1b['noise'][0, :, :]
    l1b_scaled['wavelength'] = l1b['wavelength']

    # output to netcdf file
    sim_modified_output(filename_out, l1b_scaled)

    print('=> radiometric offset added successfully ')

    return
