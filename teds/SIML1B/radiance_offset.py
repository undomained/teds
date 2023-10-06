import numpy as np
import sys
import yaml
import netCDF4 as nc
from copy import deepcopy
from ..lib.libWrite import writevariablefromname

#   main program ##############################################################


def get_l1b(path, filename):
    # getting l1b data from file
    file = path+filename+'.nc'
    input = nc.Dataset(file, mode='r')
    l1b_data = {}
    l1b_data['sza'] = deepcopy(input['GEOLOCATION_DATA']['sza'][:])
    l1b_data['saa'] = deepcopy(input['GEOLOCATION_DATA']['saa'][:])
    l1b_data['vza'] = deepcopy(input['GEOLOCATION_DATA']['vza'][:])
    l1b_data['vaa'] = deepcopy(input['GEOLOCATION_DATA']['vaa'][:])
    l1b_data['latitude']   = deepcopy(input['GEOLOCATION_DATA']['lat'][:])
    l1b_data['longitude']  = deepcopy(input['GEOLOCATION_DATA']['lon'][:])
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelengths'][:])
    l1b_data['radiance']   = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    l1b_data['noise']      = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])
    input.close()
    return l1b_data


def get_gm_data(path, filename):
    file = path+filename+'.nc'
    input = nc.Dataset(file, mode='r')
    gm_data = {}
    gm_data['sza'] = deepcopy(input['sza'][:, :])
    gm_data['saa'] = deepcopy(input['saa'][:, :])
    gm_data['vza'] = deepcopy(input['vza'][:, :])
    gm_data['vaa'] = deepcopy(input['vaa'][:, :])
    gm_data['lat'] = deepcopy(input['lat'][:, :])
    gm_data['lon'] = deepcopy(input['lon'][:, :])
    input.close()
    return gm_data


def sim_modified_output(path, filename, gm_data, l1b_output):
    nwav = l1b_output['wavelength'][0, :].size
    nact = l1b_output['radiance'][0, :, 0].size
    nalt = l1b_output['radiance'][:, 0, 0].size

    file = path+filename+'.nc'
    output = nc.Dataset(file, mode='w')
    output.title = 'Tango Carbon Level-1 B'
    output.createDimension('bins_spectral', nwav)        # spectral axis
    output.createDimension('bins_across_track', nact)    # across track axis
    output.createDimension('bins_along_track', nalt)     # along track axis

    # first geometry data
    geo_data = output.createGroup('GEOLOCATION_DATA')
    # define dimensions
    _dims = ('bins_along_track', 'bins_across_track')
    # solar zenith angle
    _ = writevariablefromname(geo_data, "solarzenithangle", _dims, gm_data['sza'])
    # solar azimuth angle
    _ = writevariablefromname(geo_data, "solarazimuthangle", _dims, gm_data['saa'])
    # viewing zenith angle
    _ = writevariablefromname(geo_data, "viewingzenithangle", _dims, gm_data['vza'])
    # viewing azimuth angle
    _ = writevariablefromname(geo_data, "viewingazimuthangle", _dims, gm_data['vaa'])
    # latitude
    _ = writevariablefromname(geo_data, "latitude", _dims, gm_data['lat'])
    # longitude
    _ = writevariablefromname(geo_data, "longitude", _dims, gm_data['lon'])

    # observation data
    obs_data = output.createGroup('OBSERVATION_DATA')
    # second radiometric data
    _dims = ('bins_across_track', 'bins_spectral')
    l1b_wave = np.zeros((nact, nwav))
    for iact in range(nact):
        l1b_wave[iact, :] = l1b_output['wavelength'][iact,:]
    writevariablefromname(obs_data, "radoff_wavelength", _dims, l1b_output['wavelength'])
    
    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    # observed Earth radiance
    writevariablefromname(obs_data, "radofff_radiance", _dims, l1b_output['radiance'])
    # observed Earth radiance noise
    writevariablefromname(obs_data, "radofff_radiance_noise", _dims, l1b_output['noise'])
    output.close()


# def sim_modified_output_old(path, filename, gm_data, l1b_output):

#     nwav = l1b_output['wavelength'][0,:].size
#     nact = l1b_output['radiance'][0, :, 0].size
#     nalt = l1b_output['radiance'][:, 0, 0].size

#     file = path+filename+'.nc'
#     output = nc.Dataset(file, mode='w')
#     output.title = 'Tango Carbon Level-1 B'
#     output.createDimension('bins_spectral', nwav)        # spectral axis
#     output.createDimension('bins_across_track', nact)    # across track axis
#     output.createDimension('bins_along_track', nalt)     # along track axis

#     # first geometry data

#     geo_data = output.createGroup('GEOLOCATION_DATA')

#     l1b_sza = geo_data.createVariable('sza', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_sza.units = 'degree'
#     l1b_sza.long_name = 'solar zenith angle'
#     l1b_sza.valid_min = 0.
#     l1b_sza.valid_max = 90.
#     l1b_sza.FillValue = -32767
#     l1b_sza[:, :] = gm_data['sza'][:, :]

#     l1b_saa = geo_data.createVariable('saa', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_saa.units = 'degree'
#     l1b_saa.long_name = 'solar azimuth angle'
#     l1b_saa.valid_min = 0.
#     l1b_saa.valid_max = +360.
#     l1b_saa.FillValue = -32767
#     l1b_saa[:, :] = gm_data['saa'][:, :]

#     l1b_vza = geo_data.createVariable('vza', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_vza.units = 'degree'
#     l1b_vza.long_name = 'viewing zenith angle'
#     l1b_vza.valid_min = 0.
#     l1b_vza.valid_max = 90.
#     l1b_vza.FillValue = -32767
#     l1b_vza[:, :] = gm_data['vza'][:, :]

#     l1b_vaa = geo_data.createVariable('vaa', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_vaa.units = 'degree'
#     l1b_vaa.long_name = 'viewing azimuth angle'
#     l1b_vaa.FillValue = -32767
#     l1b_vaa.valid_min = 0.
#     l1b_vaa.valid_max = +360.

#     l1b_vaa[:, :] = gm_data['vaa'][:, :]

#     l1b_lat = geo_data.createVariable('lat', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_lat.units = 'degree'
#     l1b_lat.long_name = 'latitude'
#     l1b_lat.FillValue = -32767
#     l1b_lat.valid_min = -90.
#     l1b_lat.valid_max = +90.
#     l1b_lat[:, :] = gm_data['lat'][:, :]

#     l1b_lon = geo_data.createVariable('lon', np.float64, ('bins_along_track', 'bins_across_track',))
#     l1b_lon.units = 'degree'
#     l1b_lon.long_name = 'longitude'
#     l1b_lon.FillValue = -32767
#     l1b_lon.valid_min = -90.
#     l1b_lon.valid_max = +90.
#     l1b_lon[:, :] = gm_data['lon'][:, :]

#     obs_data = output.createGroup('OBSERVATION_DATA')

#     # second radiometric data
#     l1b_wave = obs_data.createVariable('wavelengths', np.float64, ('bins_across_track', 'bins_spectral',))
#     l1b_wave.units = 'nm'
#     l1b_wave.long_name = 'wavelengths'
#     l1b_wave.FillValue = -32767
#     l1b_wave.valid_min = 0.
#     l1b_wave.valid_max = 80000.
#     for iact in range(nact):
#         l1b_wave[iact, :] = l1b_output['wavelength'][iact,:]

#     # l1b_sun             = output.createVariable('solar irradiance model', np.float64, ('bins_spectral',))
#     # l1b_sun.units       = 'photons/(nm m2 s)'
#     # l1b_sun.long_name   = 'solar irradiance from a model (not measured)'
#     # l1b_sun.fill_value  = 1.E-50
#     # l1b_sun.valid_min   = 0.
#     # l1b_sun.valid_max   = 1.E25
#     # l1b_sun[:]         = l1b_output['solar_irradiance'][:]

#     l1b_rad = obs_data.createVariable('radiance', np.float64, ('bins_along_track',
#                                       'bins_across_track', 'bins_spectral',))
#     l1b_rad.units = 'photons/(sr nm m2 s)'
#     l1b_rad.long_name = 'observed Earth radiance'
#     l1b_rad.FillValue = -32767
#     l1b_rad.valid_min = 0.
#     l1b_rad.valid_max = 1.E25
#     l1b_rad[:] = l1b_output['radiance'][:, :, :]

#     l1b_sig = obs_data.createVariable('radiance_noise', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
#     l1b_sig.units = 'photons/(sr nm m2 s)'
#     l1b_sig.long_name = 'relative radiance noise (1-sigma)'
#     l1b_sig.FillValue = -32767
#     l1b_sig.valid_min = 0.
#     l1b_sig.valid_max = 1.E24
#     l1b_sig[:] = l1b_output['noise'][:]

#     output.close()

#     return


def radiance_offset(paths, global_config, local_config, rad_offset):
    run_id = '_'+global_config['run_id']
    l1b_path = paths.project + paths.data_interface + paths.interface_l1b

    # paths and config parameter

    l1b_path = paths.project + paths.data_interface + paths.interface_l1b
    l1b_filename = local_config['filename']['siml1b_output'] + '_'+global_config['profile']+run_id

    l1b = get_l1b(l1b_path, l1b_filename)
    # nalt = l1b['radiance'][:, 0, 0].size
    # nact = l1b['radiance'][0, :, 0].size
    # nwav = l1b['radiance'][0, 0, :].size
    nalt, nact, nwav = l1b['radiance'].shape

    l1b_scaled = {}
    nalt_scaled = rad_offset.size
    l1b_scaled['sza'] = np.empty([nalt_scaled, nact])
    l1b_scaled['saa'] = np.empty([nalt_scaled, nact])
    l1b_scaled['vza'] = np.empty([nalt_scaled, nact])
    l1b_scaled['vaa'] = np.empty([nalt_scaled, nact])
    l1b_scaled['latitude'] = np.empty([nalt_scaled, nact])
    l1b_scaled['longitude'] = np.empty([nalt_scaled, nact])
    l1b_scaled['noise'] = np.empty([nalt_scaled, nact, nwav])
    l1b_scaled['radiance'] = np.empty([nalt_scaled, nact, nwav])
    l1b_scaled['wavelength'] = np.empty([nwav])

    for ialt in range(nalt_scaled):
        for iact in range(nact):
            offset = rad_offset[ialt]*np.max(l1b['radiance'][0, iact, :])
            l1b_scaled['radiance'][ialt, iact, :] = l1b['radiance'][0, iact, :]+offset
        l1b_scaled['sza'][ialt, :] = l1b['sza'][0, :]
        l1b_scaled['saa'][ialt, :] = l1b['saa'][0, :]
        l1b_scaled['vza'][ialt, :] = l1b['vza'][0, :]
        l1b_scaled['vaa'][ialt, :] = l1b['vaa'][0, :]
        l1b_scaled['latitude'][ialt, :] = l1b['latitude'][0, :]
        l1b_scaled['longitude'][ialt, :] = l1b['longitude'][0, :]
        l1b_scaled['noise'][ialt, :, :] = l1b['noise'][0, :, :]
    l1b_scaled['wavelength'] = l1b['wavelength']

    # get geometry data
    filename = local_config['filename']['gm_input']+'_'+global_config['profile']+run_id
    gm_path = paths.project + paths.data_interface + paths.interface_gm
    gm_data = get_gm_data(gm_path, filename)

    # output to netcdf file
    sim_modified_output(l1b_path, l1b_filename, gm_data, l1b_scaled)

    print('=> radiometric offset added successfully for runid ' +
          global_config['run_id'] + ' of profile ' + global_config['profile'])

    return
