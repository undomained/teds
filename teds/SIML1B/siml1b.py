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
from ..lib import libNumTools
from ..lib.libWrite import writevariablefromname


def sparse_isrf_convolution(isrf, mask, spectrum):
    nwav = isrf[:, 0].size
    spectrum_conv = np.empty(nwav)
    for iwav in range(nwav):
        idx = mask[iwav, :]
        spectrum_conv[iwav] = isrf[iwav, idx].dot(spectrum[idx])
    return(spectrum_conv)

def get_sgm_rad_data(filename, ialt):
    input = nc.Dataset(filename, mode='r')
    sgm_data = {}
    sgm_data['wavelength line-by-line'] = input['wavelength'][:]
    sgm_data['solar irradiance line-by-line'] = input['solar_irradiance'][:]
    sgm_data['radiance line-by-line'] = input['radiance'][ialt, :, :]
    input.close()
    return(sgm_data)


def get_gm_data(filename):
    input = nc.Dataset(filename, mode='r')
    gm_data = {}
    gm_data['sza'] = deepcopy(input['sza'][:, :])
    gm_data['saa'] = deepcopy(input['saa'][:, :])
    gm_data['vza'] = deepcopy(input['vza'][:, :])
    gm_data['vaa'] = deepcopy(input['vaa'][:, :])
    gm_data['lat'] = deepcopy(input['lat'][:, :])
    gm_data['lon'] = deepcopy(input['lon'][:, :])
    input.close()
    return (gm_data)


def sim_output(filename, gm_data, l1b_output):
    nwav = l1b_output['wavelength'].size
    nalt, nact, _ = l1b_output['radiance'].shape
    # open file
    output = nc.Dataset(filename, mode='w')
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
        l1b_wave[iact, :] = l1b_output['wavelength'][:]
    writevariablefromname(obs_data, "wavelength", _dims, l1b_output['wavelength'])

    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    # observed Earth radiance
    writevariablefromname(obs_data, "radiance", _dims, l1b_output['radiance'])
    # observed Earth radiance noise
    writevariablefromname(obs_data, "radiance_noise", _dims, l1b_output['radiance_noise'])
    # radiance error mask
    # writevariablefromname(obs_data, "radiance_mask", _dims, l1b_output['radiance_mask'])

    output.close()

#   main program ##############################################################


def simplified_instrument_model_and_l1b_processor(config):

    
    # get geometry data

    gm_data = get_gm_data(config['gm_input'])

    # target wavelengths grid
    l1b_output = {}
    l1b_output['wavelength'] = np.arange(config['spec_settings']['wave_start'],
                                         config['spec_settings']['wave_end'],
                                         config['spec_settings']['dwave'])  # nm

    # basic dimensions
    nwav = l1b_output['wavelength'].size
    nalt = gm_data['sza'][:, 0].size
    nact = gm_data['sza'][0, :].size

    # measurement array
    ymeas = np.empty([nalt, nact, nwav])

    # get line-by-line spectral grid and define some pointers

    sgm_data = get_sgm_rad_data(config['sgm_input'], ialt=0)
    wave_lbl = sgm_data['wavelength line-by-line']
    wave = l1b_output['wavelength']

    # define isrf function
    
    isrf_convolution = libNumTools.get_isrf(wave, wave_lbl, config['isrf_settings'])

    for ialt in tqdm(range(nalt)):

        # get lbl data from sgm file for scan line ialt
        sgm_data = get_sgm_rad_data(config['sgm_input'], ialt)

        for iact in range(nact):
            spectrum_lbl = np.array(sgm_data['radiance line-by-line'][iact, :])
            # isrf convolution
            ymeas[ialt, iact, :] = isrf_convolution(spectrum_lbl)

    # noise model based on SNR = a I / (sqrt (aI +b )) ; a[(e- m2 sr s nm) / phot], and [b] = e-
    snr = config['snr_model']['a_snr'] * ymeas / \
        (np.sqrt(config['snr_model']['a_snr']*ymeas + config['snr_model']['b_snr']))

    l1b_output['radiance_noise'] = ymeas/snr  # units [1]

    # random noise
    # Get nwav random numbers that are normally distributed with standard deviation 1
    np.random.seed(config['snr_model']['seed'])

    ynoise = np.empty([nalt, nact, nwav])
    for ialt in range(nalt):
        for iact in range(nact):
            noise_dis = np.random.normal(0., 1., nwav)
            # noise contribution
            ynoise[ialt, iact, :] = 1./snr[ialt, iact, :]*noise_dis*ymeas[ialt, iact, :]
    if(config['sim_with_noise']):
        l1b_output['radiance'] = ymeas+ynoise
    else:
        l1b_output['radiance'] = ymeas

    l1b_output['radiance_mask'] = np.zeros(nwav, dtype=bool)
    # output to netcdf file

    sim_output(config['l1b_output'], gm_data, l1b_output)

    print('=>siml1b calculation finished successfully')
    return
