# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# =============================================================================
#     scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

import numpy as np
import sys
import os
import pickle
import netCDF4 as nc
import yaml
from tqdm import tqdm
from copy import deepcopy
from lib.libWrite import writevariablefromname

class Dict2Class:
    """
    Convert a dictionaly to a class
    """

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


def sgm_output_radio(filename_rad, rad_output):
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
    if (gases is not None) & ("co2" in gases):
        var_co2.setncattr("source", meteodata.__getattribute__("co2_source"))
        var_co2.setncattr("emission_kgps", meteodata.__getattribute__("co2_emission_in_kgps"))
    # column_ch4
    var_ch4 = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
    if (gases is not None) & ("ch4" in gases):
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


def scene_generation_module(config):
    """
    Parameters
    ----------
    global_config : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from lib import libSGM
    from lib import libATM
    from lib import libRT
    from lib import libSURF

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
        for iscen in range(config['scene_spec']['numb_atm']+1):
            outofrange = (config['scene_spec']['scene_trans_index'][iscen] > 100) & \
                (config['scene_spec']['scene_trans_index'][iscen] < 0) 
            if (outofrange):
                sys.exit('config parameter scene_trans_index out of range')

        for iscen in range(config['scene_spec']['numb_atm']):
            ind_start = config['scene_spec']['scene_trans_index'][iscen]
            ind_end = config['scene_spec']['scene_trans_index'][iscen+1]
            albedo[0, ind_start:ind_end] = config['scene_spec']['albedo'][iscen]

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
    sgm_output_radio(config['rad_output'], rad_output)

    print('=>sgm calcultion finished successfully')
    return


if __name__ == '__main__':
    
    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )


    config = yaml.safe_load(open(sys.argv[1]))
    scene_generation_module(config)
