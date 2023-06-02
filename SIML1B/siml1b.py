# This routine provides a simplified instrument model and L1B processor, which converts the
# SGM output to a level 1B product It does not include any details on the instrument.
import numpy as np
import sys
import os
import yaml
import matplotlib.pyplot as plt
from copy import deepcopy
import netCDF4 as nc
from tqdm import tqdm

from end_to_end.lib import libNumTools


def sparse_isrf_convolution(isrf, mask, spectrum):
    nwav = isrf[:, 0].size
    spectrum_conv = np.empty(nwav)
    for iwav in range(nwav):
        idx = mask[iwav, :]
        spectrum_conv[iwav] = isrf[iwav, idx].dot(spectrum[idx])
    return(spectrum_conv)


def get_sgm_rad_data(path, filename, ialt):

    file = path+filename+'.nc'
    input = nc.Dataset(file, mode='r')

    sgm_data = {}
    sgm_data['wavelength line-by-line'] = input['wavelength'][:]
    sgm_data['solar irradiance line-by-line'] = input['solar_irradiance'][:]
    sgm_data['radiance line-by-line'] = input['radiance'][ialt, :, :]

    input.close()

    return(sgm_data)


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

    return(gm_data)

def sim_output(path, filename, gm_data, l1b_output):

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

    geo_data = output.createGroup('GEOLOCATION_DATA')

    l1b_sza = geo_data.createVariable('sza', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_sza.units = 'degree'
    l1b_sza.long_name = 'solar zenith angle'
    l1b_sza.valid_min = 0.
    l1b_sza.valid_max = 90.
    l1b_sza.FillValue = -32767
    l1b_sza[:, :] = gm_data['sza'][:, :]

    l1b_saa = geo_data.createVariable('saa', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_saa.units = 'degree'
    l1b_saa.long_name = 'solar azimuth angle'
    l1b_saa.valid_min = 0.
    l1b_saa.valid_max = +360.
    l1b_saa.FillValue = -32767
    l1b_saa[:, :] = gm_data['saa'][:, :]

    l1b_vza = geo_data.createVariable('vza', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_vza.units = 'degree'
    l1b_vza.long_name = 'viewing zenith angle'
    l1b_vza.valid_min = 0.
    l1b_vza.valid_max = 90.
    l1b_vza.FillValue = -32767
    l1b_vza[:, :] = gm_data['vza'][:, :]

    l1b_vaa = geo_data.createVariable('vaa', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_vaa.units = 'degree'
    l1b_vaa.long_name = 'viewing azimuth angle'
    l1b_vaa.FillValue = -32767
    l1b_vaa.valid_min = 0.
    l1b_vaa.valid_max = +360.
    
    l1b_vaa[:, :] = gm_data['vaa'][:, :]

    l1b_lat = geo_data.createVariable('lat', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_lat.units = 'degree'
    l1b_lat.long_name = 'latitude'
    l1b_lat.FillValue = -32767
    l1b_lat.valid_min = -90.
    l1b_lat.valid_max = +90.
    l1b_lat[:, :] = gm_data['lat'][:, :]

    l1b_lon = geo_data.createVariable('lon', np.float64, ('bins_along_track', 'bins_across_track',))
    l1b_lon.units = 'degree'
    l1b_lon.long_name = 'longitude'
    l1b_lon.FillValue = -32767
    l1b_lon.valid_min = -90.
    l1b_lon.valid_max = +90.
    l1b_lon[:, :] = gm_data['lon'][:, :]

    obs_data = output.createGroup('OBSERVATION_DATA')

    # second radiometric data
    l1b_wave = obs_data.createVariable('wavelengths', np.float64, ('bins_across_track', 'bins_spectral',))
    l1b_wave.units = 'nm'
    l1b_wave.long_name = 'wavelengths'
    l1b_wave.FillValue = -32767
    l1b_wave.valid_min = 0.
    l1b_wave.valid_max = 80000.
    for ialt in range(nalt):
        for iact in range(nact):
            l1b_wave[iact, :] = l1b_output['wavelength'][:]

    # l1b_sun             = output.createVariable('solar irradiance model', np.float64, ('bins_spectral',))
    # l1b_sun.units       = 'photons/(nm m2 s)'
    # l1b_sun.long_name   = 'solar irradiance from a model (not measured)'
    # l1b_sun.fill_value  = 1.E-50
    # l1b_sun.valid_min   = 0.
    # l1b_sun.valid_max   = 1.E25
    # l1b_sun[:]         = l1b_output['solar_irradiance'][:]

    l1b_rad = obs_data.createVariable('radiance', np.float64, ('bins_along_track',
                                      'bins_across_track', 'bins_spectral',))
    l1b_rad.units = 'photons/(sr nm m2 s)'
    l1b_rad.long_name = 'observed Earth radiance'
    l1b_rad.FillValue = -32767
    l1b_rad.valid_min = 0.
    l1b_rad.valid_max = 1.E25
    l1b_rad[:] = l1b_output['radiance'][:, :, :]

    l1b_sig = obs_data.createVariable('radiance_noise', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
    l1b_sig.units = 'photons/(sr nm m2 s)'
    l1b_sig.long_name = 'relative radiance noise (1-sigma)'
    l1b_sig.FillValue = -32767
    l1b_sig.valid_min = 0.
    l1b_sig.valid_max = 1.E24
    l1b_sig[:] = l1b_output['radiance_noise'][:]

    output.close()

    return

#   main program ##############################################################


def simplified_instrument_model_and_l1b_processor(paths, global_config, local_config):

    sgm_path = paths.project + paths.data_interface + paths.interface_sgm 

    run_id = '_'+global_config['run_id']
    
    # get geometry data
    filename = local_config['filename']['gm_input']+'_'+global_config['profile']+run_id
    gm_path = paths.project + paths.data_interface + paths.interface_gm
    gm_data = get_gm_data(gm_path, filename)

    # target wavelengths grid
    l1b_output = {}
    l1b_output['wavelength'] = np.arange(local_config['spec_settings']['wave_start'],
                                         local_config['spec_settings']['wave_end'],
                                         local_config['spec_settings']['dwave'])  # nm

    # basic dimensions
    nwav = l1b_output['wavelength'].size
    nalt = gm_data['sza'][:, 0].size
    nact = gm_data['sza'][0, :].size

    # measurement array
    ymeas = np.empty([nalt, nact, nwav])

    # get line-by-line spectral grid and define some pointers

    filename = local_config['filename']['sgm_input']+'_'+global_config['profile']+run_id
    sgm_data = get_sgm_rad_data(sgm_path, filename, ialt=0)
    wave_lbl = sgm_data['wavelength line-by-line']
    wave = l1b_output['wavelength']

    # define isrf object
    isrf = libNumTools.isrfct(wave, wave_lbl)
    isrf.get_isrf(local_config['isrf_settings'])

    for ialt in tqdm(range(nalt)):

        # get lbl data from sgm file for scan line ialt
        sgm_data = get_sgm_rad_data(sgm_path, filename, ialt)

        for iact in range(nact):
            spectrum_lbl = sgm_data['radiance line-by-line'][iact, :]
            # isrf convolution
            ymeas[ialt, iact, :] = isrf.isrf_convolution(spectrum_lbl)

    # noise model based on SNR = a I / (sqrt (aI +b )) ; a[(e- m2 sr s nm) / phot], and [b] = e-
    snr = local_config['snr_model']['a_snr'] * ymeas / \
        (np.sqrt(local_config['snr_model']['a_snr']*ymeas + local_config['snr_model']['b_snr']))

    l1b_output['radiance_noise'] = ymeas/snr  # units [1]

    # random noise
    # Get nwav random numbers that are normally distributed with standard deviation 1
    np.random.seed(local_config['snr_model']['seed'])

    ynoise = np.empty([nalt, nact, nwav])
    for ialt in range(nalt):
        for iact in range(nact):
            noise_dis = np.random.normal(0., 1., nwav)
            # noise contribution
            ynoise[ialt, iact, :] = 1./snr[ialt, iact, :]*noise_dis*ymeas[ialt, iact, :]
    if(local_config['sim_with_noise']):
        l1b_output['radiance'] = ymeas+ynoise
    else:
        l1b_output['radiance'] = ymeas

    # output to netcdf file
    siml1b_path = paths.project + paths.data_interface + paths.interface_l1b
    filename = local_config['filename']['siml1b_output'] + '_'+global_config['profile']+run_id
    print(siml1b_path)
    print(filename)
    sim_output(siml1b_path, filename, gm_data, l1b_output)

    print('=>sim calcultion finished successfully for runid ' +
          global_config['run_id'] + ' of profile ' + global_config['profile'])
    return
