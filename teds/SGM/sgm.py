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
from copy import deepcopy

import netCDF4 as nc
import numpy as np
import yaml
from scipy.interpolate import RegularGridInterpolator, griddata
from tqdm import tqdm
from lib.libWrite import writevariablefromname
from scipy.interpolate import griddata

from lib import libATM, libNumTools, libRT, libSGM, libSURF
from lib.libWrite import writevariablefromname


class Emptyclass:
    """Empty class. Data container."""
    
    pass

class Dict2Class:
    """Convert a dictionaly to a class."""

    def __init__(self, arg_dict):
        self.__dict__.update(arg_dict)


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
    return gm_data


def sgm_output_rad(filename_rad, rad_output):
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
    _ = writevariablefromname(output_rad, 'radiance_sgm', _dims, rad_output['radiance'])
    output_rad.close()


def sgm_output_atm(filename_atm, atm, albedo, gm_data, meteodata=None, gases=None):
    nalt, nact = gm_data["lat"].shape
    # write atmosphere
    nlay, nlev = atm[0][0].zlay.size, atm[0][0].zlev.size
    # file
    output_atm = nc.Dataset(filename_atm, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    # write 3d data
    # write the variables to an array
    zlev = np.zeros((nalt, nact, nlev))
    zlay = np.zeros((nalt, nact, nlay))
    dcolco2 = np.zeros((nalt, nact, nlay))
    dcolch4 = np.zeros((nalt, nact, nlay))
    dcolh2o = np.zeros((nalt, nact, nlay))
    for ialt in range(nalt):
        for iact in range(nact):
            zlay[ialt, iact, :] = atm[ialt][iact].zlay
            zlev[ialt, iact, :] = atm[ialt][iact].zlev
            dcolco2[ialt, iact, :] = atm[ialt][iact].CO2[:]
            dcolch4[ialt, iact, :] = atm[ialt][iact].CH4[:]
            dcolh2o[ialt, iact, :] = atm[ialt][iact].H2O[:]

    _dims = ('bins_along_track', 'bins_across_track', 'number_layers')
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims, zlay)
    # columndensity_co2
    _ = writevariablefromname(output_atm, 'subcol_density_co2', _dims, dcolco2)
    # columndensity_ch4
    _ = writevariablefromname(output_atm, 'subcol_density_ch4', _dims, dcolch4)
    # columndensity_h2o
    _ = writevariablefromname(output_atm, 'subcol_density_h2o', _dims, dcolh2o)
    # level height
    _dims = ('bins_along_track', 'bins_across_track', 'number_levels')
    _ = writevariablefromname(output_atm, 'levelheight', _dims, zlev)

    # Total column integrated values of CO2, CH4, H2O, and air
    xco2 = np.zeros((nalt, nact))
    xch4 = np.zeros((nalt, nact))
    xh2o = np.zeros((nalt, nact))
    col_air = np.zeros((nalt, nact))
    for ialt in range(nalt):
        for iact in range(nact):
            XAIR = np.sum(atm[ialt, iact].air[:])
            xco2[ialt, iact] = np.sum(atm[ialt, iact].CO2[:])/XAIR*1.e6  # [ppmv]
            xch4[ialt, iact] = np.sum(atm[ialt, iact].CH4[:])/XAIR*1.e9  # [ppbv]
            xh2o[ialt, iact] = np.sum(atm[ialt, iact].H2O[:])/XAIR*1.e6  # [ppmv]
            col_air[ialt, iact] = XAIR
    _dims = ('bins_along_track', 'bins_across_track')
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', _dims, albedo)
    # column_co2
    var_co2 = writevariablefromname(output_atm, 'XCO2', _dims, xco2)
    # write new attributes
    if gases is not None:
        if ("co2" in gases):
            var_co2.setncattr("source", meteodata.__getattribute__("co2_source"))
            var_co2.setncattr("emission_kgps", meteodata.__getattribute__("co2_emission_in_kgps"))
    # column_ch4
    var_ch4 = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
    if (gases is not None):
        if "ch4" in gases:
            var_ch4.setncattr("source", meteodata.__getattribute__("ch4_source"))
            var_ch4.setncattr("emission_kgps", meteodata.__getattribute__("ch4_emission_in_kgps"))

    # column_h2o
    _ = writevariablefromname(output_atm, 'XH2O', _dims, xh2o)
    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims, col_air)

    # add coordinates to SGM atmosphere
    _ = writevariablefromname(output_atm, 'latitude', _dims, gm_data["lat"])
    _ = writevariablefromname(output_atm, 'longitude', _dims, gm_data["lon"])
    output_atm.close()


def sgm_output_atm_raw(filename_atm, meteo, atm_std, gases=None):
    atm = libATM.combine_meteo_standard_atm(meteo, atm_std, gases, '_raw')
    nalt, nact = meteo.lat.shape

    # write atmosphere
    nlay, nlev = atm[0][0].zlay.size, atm[0][0].zlev.size
    # file
    output_atm = nc.Dataset(filename_atm, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    # write 3d data
    # write the variables to an array
    zlev = np.zeros((nalt, nact, nlev))
    zlay = np.zeros((nalt, nact, nlay))
    dcolco2 = np.zeros((nalt, nact, nlay))
    dcolch4 = np.zeros((nalt, nact, nlay))
    dcolh2o = np.zeros((nalt, nact, nlay))
    for ialt in range(nalt):
        for iact in range(nact):
            zlay[ialt, iact, :] = atm[ialt][iact].zlay
            zlev[ialt, iact, :] = atm[ialt][iact].zlev
            dcolco2[ialt, iact, :] = atm[ialt][iact].CO2[:]
            dcolch4[ialt, iact, :] = atm[ialt][iact].CH4[:]
            dcolh2o[ialt, iact, :] = atm[ialt][iact].H2O[:]


    _dims = ('bins_along_track', 'bins_across_track', 'number_layers')
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims, zlay)
    # columndensity_co2
    _ = writevariablefromname(output_atm, 'subcol_density_co2', _dims, dcolco2)
    # columndensity_ch4
    _ = writevariablefromname(output_atm, 'subcol_density_ch4', _dims, dcolch4)
    # columndensity_h2o
    _ = writevariablefromname(output_atm, 'subcol_density_h2o', _dims, dcolh2o)
    # level height
    _dims = ('bins_along_track', 'bins_across_track', 'number_levels')
    _ = writevariablefromname(output_atm, 'levelheight', _dims, zlev)

    # Total column integrated values of CO2, CH4, H2O, and air
    xco2 = np.zeros((nalt, nact))
    xch4 = np.zeros((nalt, nact))
    xh2o = np.zeros((nalt, nact))
    col_air = np.zeros((nalt, nact))
    for ialt in range(nalt):
        for iact in range(nact):
            XAIR = np.sum(atm[ialt, iact].air[:])
            xco2[ialt, iact] = np.sum(atm[ialt, iact].CO2[:])/XAIR*1.e6  # [ppmv]
            xch4[ialt, iact] = np.sum(atm[ialt, iact].CH4[:])/XAIR*1.e9  # [ppbv]
            xh2o[ialt, iact] = np.sum(atm[ialt, iact].H2O[:])/XAIR*1.e6  # [ppmv]
            col_air[ialt, iact] = XAIR

    albedo = meteo.albedo_raw
    _dims = ('bins_along_track', 'bins_across_track')
    
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', _dims, albedo)

    xdistance = meteo.x_new
    ydistance = meteo.y_new
    _ = writevariablefromname(output_atm, 'xdistance', 'bins_across_track', xdistance)
    _ = writevariablefromname(output_atm, 'ydistance', 'bins_along_track', ydistance)

    # column_co2
    var_co2 = writevariablefromname(output_atm, 'XCO2', _dims, xco2)
    # write new attributes
    
    if gases is not None:
        if ("co2" in gases):
            var_co2.setncattr("source", meteo.__getattribute__("co2_source"))
            var_co2.setncattr("emission_kgps", meteo.__getattribute__("co2_emission_in_kgps"))
    # column_ch4
    var_ch4 = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
    if (gases is not None):
        if "ch4" in gases:
            var_ch4.setncattr("source", meteo.__getattribute__("ch4_source"))
            var_ch4.setncattr("emission_kgps", meteo.__getattribute__("ch4_emission_in_kgps"))

    # column_h2o
    _ = writevariablefromname(output_atm, 'XH2O', _dims, xh2o)
    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims, col_air)

    # longitude/latitude coordiantes

    _ = writevariablefromname(output_atm, 'latitude', _dims, meteo.__getattribute__("lat"))
    _ = writevariablefromname(output_atm, 'longitude', _dims, meteo.__getattribute__("lon"))
    
    output_atm.close()
    return


def scene_generation_module(config):
    """
    Scent generation algorithm.
    """
    
    # first get the geometry data
    gm_data = get_gm_data(config['gm_input'])
    nact = gm_data['sza'][0].size
    nalt = len(gm_data['sza'])

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    albedo = np.zeros([nalt, nact])

    if (config['profile'] == 'individual_spectra'):
        albedo[0, :] = config['scene_spec']['albedo'][:]

    if (config['profile'] == 'single_swath'):
        ncheck = len(config['scene_spec']['albedo'])
        if (ncheck != nact):
            sys.exit("input error in sgm, nact!=100")

        albedo[0, :] = config['scene_spec']['albedo'][:]

    if (config['profile'] == 'orbit'):

        # get collocated S2 data

        file_exists = os.path.isfile(config['S2_dump'])
        if (file_exists and (not config['s2_forced'])):
            albedo = np.load(config['S2_dump'])
        else:
            print('S2 data')
            albedo = libSGM.get_sentinel2_albedo(gm_data, config)
            np.save(config['S2_dump'], albedo)    # .npy extension is added if not given
    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    nlay = config['atmosphere']['nlay']  # number of layers
    dzlay = config['atmosphere']['dzlay']
    # we assume the same standard atmosphere for all pixels of the granule
    atm_std = libATM.get_AFGL_atm_homogenous_distribution(config['afgl_input'], nlay, dzlay)

    if ((config['profile'] == 'individual_spectra') or (config['profile'] == 'single_swath')):
        xco2 = np.sum(atm_std.CO2) / np.sum(atm_std.air) * 1.E6
        atm = np.ndarray((nalt, nact), np.object_)
        for ialt in range(nalt):
            for iact in range(nact):
                atm[ialt, iact] = deepcopy(atm_std)

    if (config['profile'] == 'orbit'):
        if (config['only_afgl']):
            atm = atm_std
        else:
            # get collocated meteo data
            if ((not os.path.exists(config['meteo_dump'])) or config['meteo_forced']):
                meteodata = libATM.get_atmosphericdata(gm_data['lat'], gm_data['lon'], config['meteo'],
                                                       config['kernel_parameter'])
                # Dump meteodata dictionary into temporary pkl file
                pickle.dump(meteodata.__dict__, open(config['meteo_dump'], 'wb'))
            else:
                # Read meteodata from pickle file
                meteodata = Dict2Class(pickle.load(open(config['meteo_dump'], 'rb')))
                # combine the meteo data
            atm = libATM.combine_meteo_standard_atm(meteodata, atm_std, config["meteo"]['gases'])

    # =============================================================================
    # Write atmosphere and albedo data
    if (config['profile'] == 'orbit') & ~(config['only_afgl']):
        sgm_output_atm(config['geo_output'], atm, albedo, gm_data, meteodata, config["meteo"]['gases'])
    else:
        sgm_output_atm(config['geo_output'], atm, albedo, gm_data)

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

    nalt_ref = np.intp(nalt/2 - 0.5)
    nact_ref = np.intp(nact/2 - 0.5)

    optics = libRT.optic_abs_prop(rad_output['wavelength_lbl'], atm[nalt_ref, nact_ref].zlay)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(rad_output['wavelength_lbl'])

    # If pickle file exists read from file
    # os.path.exists(xsec_file) or conf['xsec_forced']
    if ((not os.path.exists(config['xsec_dump'])) or config['xsec_forced']):
        molec.get_data_HITRAN(config['hapi_path'], iso_ids)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm[nalt_ref][nact_ref])
        # Dump optics.prop dictionary into temporary pkl file
        pickle.dump(optics.prop, open(config['xsec_dump'], 'wb'))
    else:
        # Read optics.prop dictionary from pickle file
        optics.prop = pickle.load(open(config['xsec_dump'], 'rb'))

    # solar irradiance spectrum
    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['sun_reference'])
    rad_output['solar irradiance'] = np.interp(rad_output['wavelength_lbl'], sun['wl'], sun['phsm2nm'])

    # Calculate surface data
    surface = libSURF.surface_prop(rad_output['wavelength_lbl'])

    rad = np.empty([nalt, nact, nwav])
    for ialt in tqdm(range(nalt)):
        for iact in range(nact):
            optics.set_opt_depth_species(atm[ialt, iact], ['molec_01', 'molec_32', 'molec_07'])
            # Earth radiance spectra
            alb = [albedo[ialt, iact]]
            surface.get_albedo_poly(alb)
            rad[ialt, iact, :] = libRT.transmission(
                rad_output['solar irradiance'],
                optics,
                surface,
                np.cos(gm_data['sza'][ialt, iact]/180.*np.pi),  # mu0
                np.cos(gm_data['vza'][ialt, iact]/180.*np.pi))  # muv
    rad_output['radiance'] = rad

    # =============================================================================
    # sgm output to radiometric file
    # =============================================================================
    sgm_output_rad(config['rad_output'], rad_output)

    print('=>sgm calcultion finished successfully')
    return
 

def create_default_atmosphere(nalt, nact, atm_std):
    """Create default atmosphere.

    Parameters
    ----------
    nalt : Int
        Size across long track
    nact : Int
        Size across cross track
    atm_std : class
        Atandard atmosphere definition
    """
    atm = np.ndarray((nalt, nact), np.object_) 
    for ialt in range(nalt):                   
        for iact in range(nact):               
            atm[ialt, iact] = deepcopy(atm_std)
    return atm



def interpolate_data_irregular(meteodata, gm_data, gases):
    """Interpolate irregular data in lat-lon coordinates.

    Parameters
    ----------
    meteodata : Class 
        Class containing meterological data
    gm_data : Dict
        Dictionary containing the parameters from geometry module.
    config : List
        Contains lits of gases to be processed
    """
    # create a new container for SGM meteo data
    sgmmeteo = Emptyclass()
    sgmmeteo.__setattr__("lat", gm_data["lat"])
    sgmmeteo.__setattr__("lon", gm_data["lon"])
    sgmmeteo.__setattr__("zlev", meteodata.zlev)
    sgmmeteo.__setattr__("zlay", meteodata.zlay)

    print('Interpolating data to GM mesh...')

    dim_alt, dim_act = sgmmeteo.lat.shape   # dimensions
    dxdy = np.column_stack((meteodata.lat.ravel(), meteodata.lon.ravel()))

    # Interpolate values to GM grid
    albedo = griddata(dxdy, meteodata.albedo_conv.ravel(), (sgmmeteo.lat, sgmmeteo.lon), fill_value=0.0)
    for gas in gases:
        conv_gas = meteodata.__getattribute__(gas+"_conv")
        interpdata = np.zeros([dim_alt, dim_act, meteodata.zlay.size])
        for iz in tqdm(range(meteodata.zlay.size)):
            interpdata[:, :, iz] = griddata(dxdy, conv_gas[:,:, iz].ravel(), (sgmmeteo.lat, sgmmeteo.lon), fill_value=0.0)
        sgmmeteo.__setattr__(gas, interpdata)
    print('                     ...done')
    return albedo, sgmmeteo


def interpolate_data(meteodata, gm_data, gases):
    """Interpolate data on regular cartesian grid.

    Meteo data is in regular x-y grid and this is used to interpolate
    the values to x-y of gm grid. This is implementation is faster
    compare to interpolating data in lat-lon grid.

    Parameters
    ----------
    meteodata : Class
        Meteo data
    gm_data : Dict
        Parameters from geometry module
    gases : List
        List of gases to be processed

    Returns
    -------
    albedo : Matrix
        Albedo on gm grid
    sgmmeteo: Class
        Meteo data on gm grid
    """
    print('Interpolating data to GM mesh...')
    sgmmeteo = Emptyclass()
    sgmmeteo.__setattr__("lat", gm_data["lat"])
    sgmmeteo.__setattr__("lon", gm_data["lon"])
    sgmmeteo.__setattr__("zlev", meteodata.zlev)
    sgmmeteo.__setattr__("zlay", meteodata.zlay)
    dim_alt, dim_act = sgmmeteo.lat.shape   # dimensions
    # Interpolate albedo to GM grid
    fa = RegularGridInterpolator((meteodata.y_new, meteodata.x_new), meteodata.albedo_conv, 
                                 bounds_error=False, fill_value=0.0)
    albedo = fa((meteodata.gm_y, meteodata.gm_x))
    for gas in gases:
        interpdata = np.zeros([dim_alt, dim_act, meteodata.zlay.size])
        conv_gas = meteodata.__getattribute__(gas+"_conv")
        for iz in tqdm(range(meteodata.zlay.size)):
            fa = RegularGridInterpolator((meteodata.y_new, meteodata.x_new), conv_gas[:, :, iz], 
                                         bounds_error=False, fill_value=0.0)
            interpdata[:, :, iz] = fa((meteodata.gm_y, meteodata.gm_x))
        sgmmeteo.__setattr__(gas, interpdata)
    print('                     ...done')
    return albedo, sgmmeteo


def convolvedata(meteodata, config):
    """Convolve meteo and albedo data.

    Parameters
    ----------
    meteodata : Class
        Meteo data
    config : Dict
        Dict containing configuration parameters.

    """
    print('Convolution data...')
    conv_settings = libNumTools.getconvolutionparams(config['kernel_parameter'], meteodata.dx, meteodata.dy)
    
    # convolution of albedo   
    meteodata.__setattr__("albedo_conv", libNumTools.convolution_2d(meteodata.albedo_raw, conv_settings))
    
    for gas in config['meteo']['gases']:
        concgas = meteodata.__getattribute__(gas+"_raw")
        conv_gas = np.zeros_like(concgas)
        for iz in tqdm(range(meteodata.zlay.size)):
            # convolution
            conv_gas[:, :, iz] = libNumTools.convolution_2d(concgas[:, :, iz], conv_settings)
        meteodata.__setattr__(gas+"_conv", conv_gas)

    print('                     ...done')
    return meteodata


def scene_generation_module_new(config, sw_raw_geo_data_only=False):
    """Scene generation module.

    Parameters
    ----------
    config : Dict
       Dict containing configuration parameters.
    """
    # first get the geometry data
    gm_data = get_gm_data(config['gm_input'])
    nact = gm_data['sza'][0].size
    nalt = len(gm_data['sza'])

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================
    nlay = config['atmosphere']['nlay']  # number of layers
    dzlay = config['atmosphere']['dzlay']
    # we assume the same standard atmosphere for all pixels of the granule
    atm_std = libATM.get_AFGL_atm_homogenous_distribution(config['afgl_input'], nlay, dzlay)

    # individual spectra and single swath
    if (config['profile'] == 'individual_spectra') or (config['profile'] == 'single_swath'):
        albedo = np.zeros([nalt, nact])
        ncheck = len(config['scene_spec']['albedo'])
        if (ncheck != nact):
            sys.exit("input error in sgm, nact!=100")
        albedo[0, :] = config['scene_spec']['albedo'][:]
        atm = create_default_atmosphere(nalt, nact, atm_std)

    # Orbit
    if (config['profile'] == 'orbit'):
        if config['only_afgl']:
            atm = atm_std
            albedo = libSGM.get_sentinel2_albedo(gm_data, config)
            # functions to dump data
        else:     
            # meteorological data
            
            meteodata = libATM.get_atmosphericdata_new(gm_data['lat'], gm_data['lon'], config['meteo'])
            # get albedo on the microhh grid

            file_exists = os.path.isfile(config['S2_dump'])
            if (file_exists and (not config['s2_forced'])):
                s2_albedo = np.load(config['S2_dump'])
            else:
                s2_albedo = libSGM.get_sentinel2_albedo_new(meteodata.lat, meteodata.lon)
                np.save(config['S2_dump'], s2_albedo)    # .npy extension is added if not given

            meteodata.__setattr__("albedo_raw", s2_albedo)
            
            if(config['sw_geo_output_raw']):
                sgm_output_atm_raw(config['geo_output_raw'], meteodata, atm_std, config["meteo"]['gases'])

            #this is a temporary solution and requires later a full seperation
            #of geo-calculation and RTM
            
            if(sw_raw_geo_data_only):
                if __name__ == '__main__':
                    sys.exit()
                else:
                    return

            # convolution meteo data
            meteodata = convolvedata(meteodata, config)
            # interpolation meteo data to gm grid
            albedo, sgmmeteo = interpolate_data(meteodata, gm_data, config["meteo"]["gases"])
            # get collocated meteo data
            atm = libATM.combine_meteo_standard_atm(sgmmeteo, atm_std, config["meteo"]['gases'])

    sys.exit('jojo')
    # =============================================================================
    # Write atmosphere and albedo data
    if (config['profile'] == 'orbit') & ~(config['only_afgl']):
        sgm_output_atm(config['geo_output'], atm, albedo, gm_data, meteodata, config["meteo"]['gases'])
    else:
        sgm_output_atm(config['geo_output'], atm, albedo, gm_data)

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

    optics = libRT.optic_abs_prop(rad_output['wavelength_lbl'], atm[nalt_ref, nact_ref].zlay)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(rad_output['wavelength_lbl'])

    # If pickle file exists read from file
    # os.path.exists(xsec_file) or conf['xsec_forced']
    if ((not os.path.exists(config['xsec_dump'])) or config['xsec_forced']):
        molec.get_data_HITRAN(config['hapi_path'], iso_ids)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm[nalt_ref][nact_ref])
        # Dump optics.prop dictionary into temporary pkl file
        pickle.dump(optics.prop, open(config['xsec_dump'], 'wb'))
    else:
        # Read optics.prop dictionary from pickle file
        optics.prop = pickle.load(open(config['xsec_dump'], 'rb'))

    # solar irradiance spectrum
    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['sun_reference'])
    rad_output['solar irradiance'] = np.interp(rad_output['wavelength_lbl'], sun['wl'], sun['phsm2nm'])

    # Calculate surface data
    surface = libSURF.surface_prop(rad_output['wavelength_lbl'])

    rad = np.empty([nalt, nact, nwav])
    for ialt in tqdm(range(nalt)):
        for iact in range(nact):
            optics.set_opt_depth_species(atm[ialt, iact], ['molec_01', 'molec_32', 'molec_07'])
            # Earth radiance spectra
            alb = [albedo[ialt, iact]]
            surface.get_albedo_poly(alb)
            rad[ialt, iact, :] = libRT.transmission(
                rad_output['solar irradiance'],
                optics,
                surface,
                np.cos(gm_data['sza'][ialt, iact]/180.*np.pi),  # mu0
                np.cos(gm_data['vza'][ialt, iact]/180.*np.pi))  # muv
    rad_output['radiance'] = rad

    # =============================================================================
    # sgm output to radiometric file
    # =============================================================================
    sgm_output_rad(config['rad_output'], rad_output)

    print('=>sgm calcultion finished successfully')
    return


if __name__ == '__main__':
    
    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )


    config = yaml.safe_load(open(sys.argv[1]))
    scene_generation_module_new(config)

