# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

import os
import pickle
import sys
import netCDF4 as nc
import numpy as np
import yaml
from tqdm import tqdm

from ..lib import libNumTools, libRT, libSURF
from ..lib.libWrite import writevariablefromname
from ..lib.libNumTools import TransformCoords


class Emptyclass:
    """Empty class. Data container."""
    
    pass


class Dict2Class:
    """Convert a dictionaly to a class."""

    def __init__(self, arg_dict):
        self.__dict__.update(arg_dict)


def get_gm_data(filename):
    
    names = ['sza', 'saa', 'vza', 'vaa', 'lat', 'lon']
    
    input = nc.Dataset(filename, mode='r')
    
    gm_data = Emptyclass()

    for name in names:
        gm_data.__setattr__(name, input[name][:])
        
    input.close()
    
    return gm_data

def radsgm_output(filename_rad, rad_output):
    # write radiances
    nalt, nact, nlbl = rad_output['radiance'].shape
    # open file
    output_rad = nc.Dataset(filename_rad, mode='w')
    output_rad.title = 'Tango Carbon E2ES SGM radiometric scene'
    output_rad.createDimension('bins_spectral', nlbl)     # spectral axis
    output_rad.createDimension('bins_across_track', nact)     # across track axis
    output_rad.createDimension('bins_along_track', nalt)     # along track axis
    # wavelength
    _ = writevariablefromname(output_rad, 'wavelength', ('bins_spectral',), rad_output['wavelength_lbl'])
    # solar irradiance
    _ = writevariablefromname(output_rad, 'solarirradiance', ('bins_spectral',), rad_output['solar irradiance'])
    # radiance
    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    _ = writevariablefromname(output_rad, 'radiance', _dims, rad_output['radiance'])
    output_rad.close()

def sgm_output_atm_ref(filename, atm, albedo, gm_data, gases):

    #define dimensions
    
    dim_alt,dim_act, dim_lay = atm.zlay.shape
    dim_lev = atm.zlev.shape[2]
    
    output_atm = nc.Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', dim_alt)      # along track axis
    output_atm.createDimension('bins_across_track', dim_act)     # across track axis
    output_atm.createDimension('number_layers', dim_lay)         # layer axis
    output_atm.createDimension('number_levels', dim_lev)         # level axis

    _dims = ('bins_along_track', 'bins_across_track', 'number_layers')
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims, atm.zlay)
    # columndensity_co2
    _ = writevariablefromname(output_atm, 'subcol_density_co2', _dims, atm.CO2)
    # columndensity_ch4
    _ = writevariablefromname(output_atm, 'subcol_density_ch4', _dims, atm.CH4)
    # columndensity_h2o
    _ = writevariablefromname(output_atm, 'subcol_density_h2o', _dims, atm.H2O)
    # level height
    _dims = ('bins_along_track', 'bins_across_track', 'number_levels')
    _ = writevariablefromname(output_atm, 'levelheight', _dims, atm.zlev)

    xco2 = np.sum(atm.CO2,axis=2)/atm.air*1.e6  #[ppm]
    xch4 = np.sum(atm.CH4,axis=2)/atm.air*1.e9  #[ppb]
    xh2o = np.sum(atm.H2O,axis=2)/atm.air*1.e6  #[ppm]

    _dims = ('bins_along_track', 'bins_across_track')
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', _dims, albedo)
    # column_co2
    _ = writevariablefromname(output_atm, 'XCO2', _dims, xco2)
    # column_ch4
    _ = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
    # column_h2o
    _ = writevariablefromname(output_atm, 'XH2O', _dims, xh2o)
    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims, atm.air)

    # add coordinates to SGM atmosphere
    _ = writevariablefromname(output_atm, 'latitude', _dims, gm_data.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims, gm_data.lon)
    output_atm.close()
    
def get_geosgm_data(filename):

    input = nc.Dataset(filename, mode='r')
    
    names = ['col_air', 'dcol_ch4', 'dcol_co2', 'dcol_h2o', 'lat' , 'lon', 
             'XCH4', 'XCO2', 'XH2O', 'zlay', 'zlev','xpos', 'ypos']

    atm_data = Emptyclass()
    
    for name in names:
        atm_data.__setattr__(name, input[name][:])

    atm_data.__setattr__('albedo', input['albedo B11'][:])

    sources = [x for x in input.variables.keys() if 'source location' in x]
    emissions = [x for x in input.variables.keys() if 'emission' in x]

    for source in sources:
        attrib = source.removesuffix(' source location').lower()+'_src_zlatlon'
        atm_data.__setattr__(attrib, input[source][:])
        
    for emission in emissions:
        attrib = emission.removesuffix(' emission').lower()+'_src_kgps'
        atm_data.__setattr__(attrib, input[emission][:])
                
    input.close()
    
    return atm_data

def extract_atm(atm,ialt,iact):

    atm_ext = Emptyclass()

    atm_ext.__setattr__('lat',atm.lat[ialt,iact])
    atm_ext.__setattr__('lon',atm.lon[ialt,iact])
    atm_ext.__setattr__('CO2',atm.CO2[ialt,iact,:])
    atm_ext.__setattr__('CH4',atm.CH4[ialt,iact,:])
    atm_ext.__setattr__('H2O',atm.H2O[ialt,iact,:])
    atm_ext.__setattr__('zlay',atm.zlay[ialt,iact,:])
    atm_ext.__setattr__('zlev',atm.zlev[ialt,iact,:])
    atm_ext.__setattr__('air',atm.air[ialt,iact])

    return atm_ext

def Carbon_radiation_scene_generation(config):
    """Scene generation module.

    Parameters
    ----------
    config : Dict
       Dict containing configuration parameters.
    """
    #get the geometry data
    gm_data = get_gm_data(config['gm_input'])
    nalt, nact = gm_data.sza.shape
    
    #get data 
    atm_org = get_geosgm_data(config['geo_output'])
    
    #convolution with instrument spatial response
    atm_conv = libNumTools.convolvedata(atm_org, config)
 
    # create a transform method 
    trans = TransformCoords(atm_org.co2_src_zlatlon[1:])

    # convert lat-lon of gm to x-y and get bounds
    gm_data.xpos, gm_data.ypos = trans.latlon2xymts(gm_data.lat, gm_data.lon)

    #interpolate convolved data to gm grid
    albedo, atm = libNumTools.interpolate_data_regular(atm_conv, gm_data, config["conv_gases"])

    #store sgm data for which RTM is done
    sgm_output_atm_ref(config['geo_output_ref'], atm, albedo, gm_data, config["conv_gases"])

    # =============================================================================
    #  Radiative transfer simulations
    # =============================================================================

    print('radiative transfer simuation...')
    # define line-by-line wavelength grid
    rad_output = {}
    rad_output['wavelength_lbl'] = np.arange(config['spec_settings']['wave_start'], config['spec_settings']['wave_end'],
                                             config['spec_settings']['dwave'])  # nm

    nwav = len(rad_output['wavelength_lbl'])
    # generate optics object for one representative model atmosphere of the domain

    nalt_ref = np.int0(nalt/2 - 0.5)
    nact_ref = np.int0(nact/2 - 0.5)
    atm_ref  = extract_atm(atm,nalt_ref,nact_ref)
    
    optics = libRT.optic_abs_prop(rad_output['wavelength_lbl'], atm_ref.zlay)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(rad_output['wavelength_lbl'])

    # If pickle file exists read from file
    # os.path.exists(xsec_file) or conf['xsec_forced']
    if ((not os.path.exists(config['xsec_dump'])) or config['xsec_forced']):
        molec.get_data_HITRAN(config['hapi_path'], iso_ids)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm_ref)
        # Dump optics.prop dictionary into temporary pkl file
        pickle.dump(optics.prop, open(config['xsec_dump'], 'wb'))
    else:
        # Read optics.prop dictionary from pickle file
        optics.prop = pickle.load(open(config['xsec_dump'], 'rb'))

    #filename = './xsection.nc'
    #output_optics(filename, atm_std, optics)
    #sys.exit()
    
    # solar irradiance spectrum
    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['sun_reference'])
    rad_output['solar irradiance'] = np.interp(rad_output['wavelength_lbl'], sun['wl'], sun['phsm2nm'])

    # Calculate surface data
    surface = libSURF.surface_prop(rad_output['wavelength_lbl'])

    rad = np.empty([nalt, nact, nwav])
    for ialt in tqdm(range(nalt)):
        for iact in range(nact):
            atm_ext = extract_atm(atm,ialt,iact)
            optics.set_opt_depth_species(atm_ext, ['molec_01', 'molec_32', 'molec_07'])
            # Earth radiance spectra
            alb = [albedo[ialt, iact]]
            mu_sza = np.cos(np.deg2rad(gm_data.sza[ialt, iact]))
            mu_vza = np.cos(np.deg2rad(gm_data.vza[ialt, iact]))
            surface.get_albedo_poly(alb)
            rad[ialt, iact, :] = libRT.transmission(
                rad_output['solar irradiance'], optics, surface, mu_sza, mu_vza)  
            
                       
    rad_output['radiance'] = rad

    # =============================================================================
    # sgm output to radiometric file
    # =============================================================================

    radsgm_output(config['rad_output'], rad_output)

    print('=>Carbon radsgm calculation finished successfully')
    return

if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    Carbon_radiation_scene_generation(config)
    