import numpy as np
import sys
import yaml
import netCDF4 as nc

from modules.lib import libINV

#   main program ##############################################################


def sim_modified_output(path, filename, l1b_output):

    nwav = l1b_output['wavelength'].size
    nact = l1b_output['radiance'][0, :, 0].size
    nalt = l1b_output['radiance'][:, 0, 0].size

    file = path+filename+'.nc'
    output = nc.Dataset(file, mode='w')
    output.title = 'Tango Carbon Level-1 B'
    output.createDimension('bins_spectral', nwav)        # spectral axis
    output.createDimension('bins_across_track', nact)    # across track axis
    output.createDimension('bins_along_track', nalt)     # along track axis

    # first geometry data

    l1b_sza = output.createVariable('sza', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_sza.units = 'degree'
    l1b_sza.long_name = 'solar zenith angle'
    l1b_sza.fill_value = 1.E-50
    l1b_sza.valid_min = 0.
    l1b_sza.valid_max = 90.
    l1b_sza[:, :] = l1b_output['sza'][:, :]

    l1b_saa = output.createVariable('saa', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_saa.units = 'degree'
    l1b_saa.long_name = 'solar azimuth angle'
    l1b_saa.fill_value = 1.E-50
    l1b_saa.valid_min = 0.
    l1b_saa.valid_max = +360.
    l1b_saa[:, :] = l1b_output['saa'][:, :]

    l1b_vza = output.createVariable('vza', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_vza.units = 'degree'
    l1b_vza.long_name = 'viewing zenith angle'
    l1b_vza.fill_value = 1.E-50
    l1b_vza.valid_min = 0.
    l1b_vza.valid_max = 90.
    l1b_vza[:, :] = l1b_output['vza'][:, :]

    l1b_vaa = output.createVariable('vaa', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_vaa.units = 'degree'
    l1b_vaa.long_name = 'viewing azimuth angle'
    l1b_vaa.fill_value = 1.E-50
    l1b_vaa.valid_min = 0.
    l1b_vaa.valid_max = +360.
    l1b_vaa[:, :] = l1b_output['vaa'][:, :]

    l1b_lat = output.createVariable('lat', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_lat.units = 'degree'
    l1b_lat.long_name = 'latitude'
    l1b_lat.fill_value = 1.E-50
    l1b_lat.valid_min = -90.
    l1b_lat.valid_max = +90.
    l1b_lat[:, :] = l1b_output['latitude'][:]

    l1b_lon = output.createVariable('lon', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_lon.units = 'degree'
    l1b_lon.long_name = 'longitude'
    l1b_lon.fill_value = 1.E-50
    l1b_lon.valid_min = -90.
    l1b_lon.valid_max = +90.
    l1b_lon[:, :] = l1b_output['longitude'][:]

    # second radiometric data
    l1b_wave = output.createVariable('wavelength', np.float64, ('bins_spectral',))
    l1b_wave.units = 'nm'
    l1b_wave.long_name = 'wavelength'
    l1b_wave.fill_value = 1.E-50
    l1b_wave.valid_min = 0.
    l1b_wave.valid_max = 80000.
    l1b_wave[:] = l1b_output['wavelength'][:]

    # l1b_sun             = output.createVariable('solar irradiance model', np.float64, ('bins_spectral',))
    # l1b_sun.units       = 'photons/(nm m2 s)'
    # l1b_sun.long_name   = 'solar irradiance from a model (not measured)'
    # l1b_sun.fill_value  = 1.E-50
    # l1b_sun.valid_min   = 0.
    # l1b_sun.valid_max   = 1.E25
    # l1b_sun[:]         = l1b_output['solar_irradiance'][:]

    l1b_rad = output.createVariable('radiance', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
    l1b_rad.units = 'photons/(sr nm m2 s)'
    l1b_rad.long_name = 'observed Earth radiance'
    l1b_rad.fill_value = 1.E-50
    l1b_rad.valid_min = 0.
    l1b_rad.valid_max = 1.E25
    l1b_rad[:] = l1b_output['radiance'][:]

    l1b_sig = output.createVariable('noise', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
    l1b_sig.units = 'photons/(sr nm m2 s)'
    l1b_sig.long_name = 'relative radiance noise (1-sigma)'
    l1b_sig.fill_value = 1.E-50
    l1b_sig.valid_min = 0.
    l1b_sig.valid_max = 1.E24
    l1b_sig[:] = l1b_output['noise'][:]

    output.close()

    return


def radiance_offset(global_config, rad_offset):

    siml1b_path = global_config['path']['e2es_project'] + global_config['path']['siml1b_path']
    l1b_path = global_config['path']['e2es_project'] + global_config['path']['interface_data']
    config = yaml.safe_load(open(siml1b_path+'siml1b_config.yaml'))
    run_id = '_'+global_config['run_id']

    # output to netcdf file
    filename = config['path']['siml1b_output']+config['filename']['siml1b_output'] \
        + '_'+global_config['profile']+'_exp0001'

    # paths and config parameter

    l1b = libINV.get_l1b(l1b_path, filename)

    nalt = l1b['radiance'][:, 0, 0].size
    nact = l1b['radiance'][0, :, 0].size
    nwav = l1b['radiance'][0, 0, :].size

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

    # output to netcdf file
    filename = config['path']['siml1b_output']+config['filename']['siml1b_output'] \
        + '_'+global_config['profile']+run_id
    sim_modified_output(l1b_path, filename, l1b_scaled)

    print('=> radiometric offset added successfully for runid ' +
          global_config['run_id'] + ' of profile ' + global_config['profile'])

    return
