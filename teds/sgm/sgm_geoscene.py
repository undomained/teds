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
    """Read a list of Sentinel 2 albedos from a NetCDF file."""
    nc = Dataset(filename)
    albedos = []

    for group in [x for x in nc.groups if x != 'SCL']:

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


def get_sentinel2_scl(filename: str) -> DataArray:
    """Read Sentinel 2 surface classification layer from NetCDF file."""
    nc = Dataset(filename)
    grp = nc['SCL']
    scl = DataArray(grp['scl'][:],
                    dims=('y', 'x'),
                    coords={'y': grp['y'][:], 'x': grp['x'][:]})
    scl.attrs['gsd'] = grp['gsd'][:]
    scl.rio.write_crs(grp.crs, inplace=True)
    return scl


class Emptyclass:
    """Empty class. Data container."""

    pass


def get_gm_data(filename):

    input = nc.Dataset(filename, mode='r')
    gm_data = Emptyclass()
    gm_data.__setattr__('sza', input['solar_zenith'][:])
    gm_data.__setattr__('saa', input['solar_azimuth'][:])
    gm_data.__setattr__('vza', input['sensor_zenith'][:])
    gm_data.__setattr__('vaa', input['sensor_azimuth'][:])
    gm_data.__setattr__('lat', input['latitude'][:])
    gm_data.__setattr__('lon', input['longitude'][:])
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
        if gas != 'air':
            # column mixing ratio
            _ = writevariablefromname(output_atm, 'X'+gas, _dims2d,
                                      constants.__getattribute__('scale_X'+gas)*
                                      atm.__getattribute__('X'+gas))
        
#    if(config['profile']=='orbit'):

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

    var = output_atm.createVariable('latitude', 'f8', _dims2d)
    var.long_name = 'latitudes'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = atm.lat

    var = output_atm.createVariable('longitude', 'f8', _dims2d)
    var.long_name = 'longitudes'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = atm.lon

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

def geosgm_output_ind_spec(filename, atm):
    # write geophysical scene data to output

    nalt, nact, nlay = atm.zlay.shape
    nlev = nlay+1

    output_atm = nc.Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

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
    _ = writevariablefromname(output_atm, 'albedo_B11',_dims2d,atm.albedo)

    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims2d, atm.col_air)
    # longitude/latitude coordiantes
    _ = writevariablefromname(output_atm, 'latitude', _dims2d, atm.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims2d, atm.lon)
    # xpos and ypos
    _ = writevariablefromname(output_atm, 'xpos', 'bins_across_track', atm.xpos)
    _ = writevariablefromname(output_atm, 'ypos', 'bins_along_track', atm.ypos)

    output_atm.close()

    return


def geoscene_generation(config: dict) -> None:
    """Generate a geophysical scene.

    Args:
      config
        Configuration dictionary

    """
    # first  the geometry data
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
        if (len(config['scene_spec']['albedo']) != nact):
            sys.exit("input error in sgm, albedo dimension not consistent with gm")

        alb = np.zeros([nalt, nact])
        alb[0, :] = config['scene_spec']['albedo'][:]
        albedo = DataArray(alb,
                           dims=('y', 'x'),)
        albedo.attrs['gsd'] = 300  #m 
        albedo.attrs['band_label'] = 'B11'
        albedo.attrs['central_wavelength'] = 1620
        albedo.attrs['bandwidth'] = '60'

        atm = libATM.create_atmosphere_ind_spectra(nalt, nact, atm_std, albedo, gm_data)

        geosgm_output(config['io_files']['output_geo'], atm)
    # Orbit
    if (config['profile'] == 'orbit'):

        # meteorological data

        meteodata = libATM.get_atmosphericdata_new(gm_data.lat, gm_data.lon, config['io_files']['meteo'])

        # Get albedo on the microHH grid
        s2_albedos = get_sentinel2_albedo(config['io_files']['input_s2'])

        #replace nan with closest non-nan value
        for s2_alb in s2_albedos:
            mask = np.isnan(s2_alb.values)
            idx = np.where(~mask,np.arange(mask.shape[1]),0)
            np.maximum.accumulate(idx,axis=1, out=idx)
            s2_alb.values[mask] = s2_alb.values[np.nonzero(mask)[0], idx[mask]]

        s2_albedos = libSGM.interp_sentinel2_albedo(
            s2_albedos,
            meteodata.lat,
            meteodata.lon,
            config['sentinel2']['band_label'])

        meteodata.__setattr__("albedo", s2_albedos)

        atm = libATM.combine_meteo_standard_atm_new(meteodata, atm_std, config)

        geosgm_output(config['io_files']['output_geo'], atm)

    print('=>sgm geoscene calculation finished successfully')

if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    geoscene_generation(config)
