# This routine provides a simplified instrument model and L1B processor, which converts the
# SGM output to a level 1B product It does not include any details on the instrument.
import numpy as np
import sys
import os
import yaml
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
    output = nc.Dataset(filename, mode='w')
    output.title = 'Tango Carbon level 1B data'

    nalt, nact, nwav = l1b_output['radiance'].shape
    output.createDimension('wavelength', nwav)
    output.createDimension('across_track_sample', nact)
    output.createDimension('along_track_sample', nalt)

    nc_grp = output.createGroup('geolocation_data')
    _dims = ('along_track_sample', 'across_track_sample')

    nc_var = nc_grp.createVariable(
        'latitude', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'latitude at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees_north'
    nc_var[:] = gm_data['lat']

    nc_var = nc_grp.createVariable(
        'longitude', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'longitude at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees_east'
    nc_var[:] = gm_data['lon']

    nc_var = nc_grp.createVariable(
        'solar_zenith', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'solar zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = gm_data['sza']

    nc_var = nc_grp.createVariable(
        'solar_azimuth', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'solar azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = gm_data['saa']

    nc_var = nc_grp.createVariable(
        'sensor_zenith', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'sensor zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = gm_data['vza']

    nc_var = nc_grp.createVariable(
        'sensor_azimuth', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'sensor azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = gm_data['vaa']

    nc_grp = output.createGroup('observation_data')
    _dims = ('across_track_sample', 'wavelength')

    l1b_wave = np.zeros((nact, nwav))
    for iact in range(nact):
        l1b_wave[iact, :] = l1b_output['wavelength'][:]
    writevariablefromname(
        nc_grp, 'wavelength', _dims, l1b_output['wavelength'])

    _dims = ('along_track_sample', 'across_track_sample', 'wavelength')

    writevariablefromname(nc_grp, 'radiance', _dims, l1b_output['radiance'])

    nc_var = nc_grp.createVariable(
        'radiance_stdev', 'f8', _dims, fill_value=-32767.0)
    nc_var.long_name = 'relative radiance noise (1-sigma)'
    nc_var.units = 'photons/(sr nm m2 s)'
    nc_var.valid_min = 0.0
    nc_var.valid_max = 1e24
    nc_var[:] = l1b_output['radiance_noise']

    output.close()


#   main program ##############################################################
def simplified_instrument_model_and_l1b_processor(config):

    # get geometry data

    gm_data = get_gm_data(config['io_files']['input_gm'])

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

    sgm_data = get_sgm_rad_data(config['io_files']['input_sgm'], ialt=0)
    wave_lbl = sgm_data['wavelength line-by-line']
    wave = l1b_output['wavelength']

    # define isrf function

    isrf_convolution = libNumTools.get_isrf(wave, wave_lbl, config['isrf_settings'])

    for ialt in tqdm(range(nalt)):

        # get lbl data from sgm file for scan line ialt
        sgm_data = get_sgm_rad_data(config['io_files']['input_sgm'], ialt)

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
        l1b_output['radiance'] = ymeas + ynoise
    else:
        l1b_output['radiance'] = ymeas

    l1b_output['radiance_mask'] = np.zeros(nwav, dtype=bool)
    # output to netcdf file

    sim_output(config['io_files']['output_l1b'], gm_data, l1b_output)

    print('=>siml1b calculation finished successfully')
    return
