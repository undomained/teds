# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     geophysical scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

import sys
import netCDF4 as nc
import numpy as np
import yaml
from ..lib import constants
from ..lib import libATM, libSGM
from ..lib.libWrite import writevariablefromname

from netCDF4 import Dataset
from xarray import DataArray
from typing import List


def get_sentinel2_albedo(filename: str) -> List[DataArray]:
    '''Read a list of Sentinel 2 albedos from a NetCDF file.

    '''
    nc = Dataset(filename)
    albedos = []
    for group in nc.groups:
        albedo = DataArray(nc[group]['albedo'][:],
                           dims=('y', 'x'),
                           coords={
                               'y': nc[group]['y'][:],
                               'x': nc[group]['x'][:]
                           })
        albedo.attrs['gsd'] = nc[group]['gsd'][:]
        albedo.attrs['band_label'] = group
        albedo.rio.write_crs(nc[group].crs, inplace=True)
        albedos.append(albedo)
    return albedos


class Emptyclass:
    """Empty class. Data container."""

    pass


def get_gm_data(filename):

    names = ['sza', 'saa', 'vza', 'vaa', 'lat', 'lon']

    input = nc.Dataset(filename, mode='r')

    gm_data = Emptyclass()

    for name in names:
        gm_data.__setattr__(name, input[name][:])

    input.close()

    return gm_data


def geosgm_output(filename, atm):
    # write geophysical scene data to output

    nalt, nact, nlay = atm.zlay.shape
    nlev = nlay+1

    output_atm = nc.Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis
    output_atm.createDimension('location', 3)                 # source coordiantes
    output_atm.createDimension('emission', 1)                 # emission strength

    _dims3dlay = ('bins_along_track', 'bins_across_track', 'number_layers')
    _dims3dlev = ('bins_along_track', 'bins_across_track', 'number_levels')
    _dims2d    = ('bins_along_track', 'bins_across_track')

    gases = [x.removeprefix('dcol_') for x in atm.__dict__.keys() if 'dcol_' in x]

    # level height
    _ = writevariablefromname(output_atm, 'levelheight', _dims3dlev, atm.zlev)
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims3dlay, atm.zlay)

    for gas in gases:
        # subcolumn density
        _ = writevariablefromname(output_atm, 'subcol_density_'+gas,
                                  _dims3dlay, atm.__getattribute__('dcol_'+gas))
        # column mixing ratio
        _ = writevariablefromname(output_atm, 'X'+gas, _dims2d, 
                                  constants.__getattribute__('scale_X'+gas)* 
                                  atm.__getattribute__('X'+gas))

    # albedo
    for s2_albedo in atm.albedo:
        var_alb = writevariablefromname(
            output_atm,
            'albedo_' + s2_albedo.band_label,
            _dims2d,
            s2_albedo.values)
        var_alb.setncattr("central wavelength", s2_albedo.central_wavelength)
        var_alb.setncattr("band width", s2_albedo.bandwidth)

#    _ = writevariablefromname(output_atm, 'albedo', _dims2d, atm.albedo)
    # xpos and ypos       
    _ = writevariablefromname(output_atm, 'xpos', 'bins_across_track', atm.xpos)
    _ = writevariablefromname(output_atm, 'ypos', 'bins_along_track', atm.ypos)

    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims2d, atm.col_air)
    # longitude/latitude coordiantes
    _ = writevariablefromname(output_atm, 'latitude', _dims2d, atm.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims2d, atm.lon)

    #information on emission source    
    substr = 'source'
    attr = atm.__dict__.keys()
    attr_src = [string for string in attr if substr in string]
    for src in attr_src:
        _ = writevariablefromname(output_atm, src+'_location', 'location', atm.__getattribute__(src))

    substr = 'emission'
    attr = atm.__dict__.keys()
    attr_emi = [string for string in attr if substr in string]
    for emi in attr_emi:
        _ = writevariablefromname(output_atm, emi.removesuffix('_in_kgps'), 'emission' , atm.__getattribute__(emi))

    output_atm.close()
            
    return

def geoscene_generation(config):
    """Scene generation module.

    Parameters
    ----------
    config : Dict
       Dict containing configuration parameters.
    """
    # first get the geometry data
    gm_data = get_gm_data(config['io_files']['input_gm'])
    nalt, nact = gm_data.sza.shape

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================
    nlay = config['atmosphere']['nlay']  # number of layers
    dzlay = config['atmosphere']['dzlay']
    # we assume the same standard atmosphere for all pixels of the granule

    atm_std = libATM.get_AFGL_atm_homogenous_distribution(
        config['io_files']['input_afgl'], nlay, dzlay, config['scale_gas']['xco2'],
        config['scale_gas']['xch4'], config['scale_gas']['xh2o'])

    # individual spectra and single swath
    if (config['profile'] == 'individual_spectra'):
        #use this profile whe you want to study the retrieval dependence as a
        #function of one parameter
        if(nalt!= 1):
            sys.exit("input error in sgm, for profile = indiudual spectra, nalt!=1")
        albedo = np.zeros([nalt, nact])
        if (len(config['scene_spec']['albedo']) != nact):
            sys.exit("input error in sgm, albedo dimension not consistent with gm")
        albedo[0, :] = config['scene_spec']['albedo'][:]
        atm = libATM.create_atmosphere_ind_spectra(nalt, nact, atm_std, albedo)

    # Orbit
    if (config['profile'] == 'orbit'):

        # meteorological data
        meteodata = libATM.get_atmosphericdata_new(
            gm_data.lat, gm_data.lon, config['meteo'])

        meteodata = libATM.get_atmosphericdata_new(gm_data.lat, gm_data.lon, config['io_files']['meteo'])

        # Get albedo on the microHH grid
        s2_albedos = get_sentinel2_albedo(config['io_files']['input_s2'])
        
        s2_albedos = libSGM.interp_sentinel2_albedo(
            s2_albedos,
            meteodata.lat,
            meteodata.lon,
            config['kernel_parameter'],
            config['sentinel2']['band_section'])
        
        meteodata.__setattr__("albedo", s2_albedos)
        
        # get albedo on the microhh grid
#        file_exists = os.path.isfile(config['io_files']['input_s2'])
#        
#        print('file exists: ', file_exists)
#        sys.exit()
#        if (file_exists and (not config['S2']['force_data_ret'])):
#            albedo = pickle.load(open(config['S2_dump'], 'rb'))
#        else:
#            albedo = libSGM.get_sentinel2_albedo_new(meteodata.lat, meteodata.lon, config['S2']['band_section'])
#            pickle.dump(albedo, open(config['S2_dump'], 'wb')) # .pkl extension is added if not given
#            
#        meteodata.__setattr__("albedo", albedo)
        
        atm = libATM.combine_meteo_standard_atm_new(meteodata, atm_std, config)

        geosgm_output(config['io_files']['output_geo'], atm)
        
    print('=>sgm geoscene calculation finished successfully')


if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    geoscene_generation(config)
