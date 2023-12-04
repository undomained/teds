import numpy as np
import sys
import yaml
import netCDF4 as nc
from copy import deepcopy
from ..lib.libWrite import writevariablefromname

#   main program ##############################################################


def get_l1b(filename):
    # getting l1b data from file
    input = nc.Dataset(filename, mode='r')
    l1b_data = {}
    l1b_data['sza'] = deepcopy(input['GEOLOCATION_DATA']['sza'][:])
    l1b_data['saa'] = deepcopy(input['GEOLOCATION_DATA']['saa'][:])
    l1b_data['vza'] = deepcopy(input['GEOLOCATION_DATA']['vza'][:])
    l1b_data['vaa'] = deepcopy(input['GEOLOCATION_DATA']['vaa'][:])
    l1b_data['lat']   = deepcopy(input['GEOLOCATION_DATA']['lat'][:])
    l1b_data['lon']  = deepcopy(input['GEOLOCATION_DATA']['lon'][:])
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelength'][:])
    l1b_data['radiance']   = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    l1b_data['noise']      = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])
    input.close()
    return l1b_data

def sim_modified_output(filename, l1b_output):
    nwav = l1b_output['wavelength'][0, :].size
    nact = l1b_output['radiance'][0, :, 0].size
    nalt = l1b_output['radiance'][:, 0, 0].size

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
    _ = writevariablefromname(geo_data, "solarzenithangle", _dims, l1b_output['sza'])
    # solar azimuth angle
    _ = writevariablefromname(geo_data, "solarazimuthangle", _dims, l1b_output['saa'])
    # viewing zenith angle
    _ = writevariablefromname(geo_data, "viewingzenithangle", _dims, l1b_output['vza'])
    # viewing azimuth angle
    _ = writevariablefromname(geo_data, "viewingazimuthangle", _dims, l1b_output['vaa'])
    # latitude
    _ = writevariablefromname(geo_data, "latitude", _dims, l1b_output['lat'])
    # longitude
    _ = writevariablefromname(geo_data, "longitude", _dims, l1b_output['lon'])

    # observation data
    obs_data = output.createGroup('OBSERVATION_DATA')
    # second radiometric data
    _dims = ('bins_across_track', 'bins_spectral')
    l1b_wave = np.zeros((nact, nwav))
    for iact in range(nact):
        l1b_wave[iact, :] = l1b_output['wavelength'][iact,:]
    writevariablefromname(obs_data, "wavelength", _dims, l1b_output['wavelength'])
    
    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    # observed Earth radiance
    writevariablefromname(obs_data, "radiance", _dims, l1b_output['radiance'])
    # observed Earth radiance noise
    writevariablefromname(obs_data, "radiance_noise", _dims, l1b_output['noise'])
    output.close()

def radiance_offset(filename_in, filename_out, rad_offset):

    l1b = get_l1b(filename_in)
    nalt, nact, nwav = l1b['radiance'].shape

    l1b_scaled = {}
    l1b_scaled['sza'] = np.empty([nalt, nact])
    l1b_scaled['saa'] = np.empty([nalt, nact])
    l1b_scaled['vza'] = np.empty([nalt, nact])
    l1b_scaled['vaa'] = np.empty([nalt, nact])
    l1b_scaled['lat'] = np.empty([nalt, nact])
    l1b_scaled['lon'] = np.empty([nalt, nact])
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
        l1b_scaled['lat'][ialt, :] = l1b['lat'][0, :]
        l1b_scaled['lon'][ialt, :] = l1b['lon'][0, :]
        l1b_scaled['noise'][ialt, :, :] = l1b['noise'][0, :, :]
    l1b_scaled['wavelength'] = l1b['wavelength']

    # output to netcdf file
    sim_modified_output(filename_out, l1b_scaled)

    print('=> radiometric offset added successfully ')

    return
