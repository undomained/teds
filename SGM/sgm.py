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
import matplotlib.pyplot as plt
import netCDF4 as nc
import yaml
from tqdm import tqdm
from copy import deepcopy


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

    return(gm_data)


def sgm_output(filename_rad, filename_atm, rad_output, atm, albedo):

    nalt = len(rad_output['radiance'][:, 0, 0])
    nact = len(rad_output['radiance'][0, :, 0])
    nlbl = len(rad_output['radiance'][0, 0, :])

    output_rad = nc.Dataset(filename_rad, mode='w')

    output_rad.title = 'Tango Carbon E2ES SGM radiometric scene'
    output_rad.createDimension('bins_spectral', nlbl)     # spectral axis
    output_rad.createDimension('bins_across_track', nact)     # across track axis
    output_rad.createDimension('bins_along_track', nalt)     # along track axis

    sgm_wavelbl = output_rad.createVariable('wavelength', np.float64, ('bins_spectral',))
    sgm_wavelbl.units = 'nm'
    sgm_wavelbl.long_name = 'wavelength line-by-line'
    sgm_wavelbl.valid_min = 0.
    sgm_wavelbl.valid_max = 8000.
    sgm_wavelbl.FillValue = -32767
    sgm_wavelbl[:] = rad_output['wavelength lbl'][:]

    sgm_sunlbl = output_rad.createVariable('solar_irradiance', np.float64, ('bins_spectral',))
    sgm_sunlbl.units = 'photons / (nm m2 s)'
    sgm_sunlbl.long_name = 'solar irradiance line-by-line'
    sgm_sunlbl.valid_min = 0.
    sgm_sunlbl.valid_max = 1.E30
    sgm_sunlbl.FillValue = -32767
    sgm_sunlbl[:] = rad_output['solar irradiance'][:]

    sgm_radlbl = output_rad.createVariable(
        'radiance', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
    sgm_radlbl.units = 'photons / (sr nm m2 s)'
    sgm_radlbl.long_name = 'radiance line-by-line'
    sgm_radlbl.valid_min = 0.
    sgm_radlbl.valid_max = 1.E30
    sgm_radlbl.FillValue = -32767
    sgm_radlbl[:] = rad_output['radiance'][:]

    output_rad.close()

    nlay = atm[0][0].zlay.size
    nlev = atm[0][0].zlev.size

    output_atm = nc.Dataset(filename_atm, mode='w')

    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    sgm_albedo = output_atm.createVariable('albedo', np.float64, ('bins_along_track', 'bins_across_track',))
    sgm_albedo.units = '-'
    sgm_albedo.long_name = 'Lambertian surface albedo'
    sgm_albedo.valid_min = 0.
    sgm_albedo.valid_max = 1.
    sgm_albedo.FillValue = -32767
    sgm_albedo[:] = albedo[:]

    sgm_zlay = output_atm.createVariable('zlay', np.float64, ('bins_along_track',
                                                              'bins_across_track', 'number_layers',))
    sgm_zlay.units = 'm'
    sgm_zlay.long_name = 'central layer height'
    sgm_zlay.valid_min = 0.
    sgm_zlay.valid_max = 1.E+5
    sgm_zlay.FillValue = -32767
    
    for ialt in range(nalt):
        for iact in range(nact):
            sgm_zlay[ialt, iact, :] = atm[ialt][iact].zlay

    sgm_zlev = output_atm.createVariable('zlev', np.float64, ('bins_along_track',
                                                              'bins_across_track', 'number_levels',))
    sgm_zlev.units = 'm'
    sgm_zlev.long_name = 'level height'
    sgm_zlev.valid_min = 0.
    sgm_zlev.valid_max = 1.E+5
    sgm_zlev.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            sgm_zlev[ialt, iact, :] = atm[ialt][iact].zlev

    sgm_dcol_co2 = output_atm.createVariable(
        'dcol_co2', np.float64, ('bins_along_track', 'bins_across_track', 'number_layers',))
    sgm_dcol_co2.units = 'molec./cm2'
    sgm_dcol_co2.long_name = 'CO2 layer column density'
    sgm_dcol_co2.valid_max = 1.E+23
    sgm_dcol_co2.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            sgm_dcol_co2[ialt, iact, :] = atm[ialt][iact].CO2[:]

    sgm_dcol_ch4 = output_atm.createVariable(
        'dcol_ch4', np.float64, ('bins_along_track', 'bins_across_track', 'number_layers',))
    sgm_dcol_ch4.units = 'molec./cm2'
    sgm_dcol_ch4.long_name = 'CH4 layer column density'
    sgm_dcol_ch4.valid_min = 0.
    sgm_dcol_ch4.valid_max = 1.E+23
    sgm_dcol_ch4.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            sgm_dcol_ch4[ialt, iact, :] = atm[ialt][iact].CH4[:]

    sgm_dcol_h2o = output_atm.createVariable(
        'dcol_h2o', np.float64, ('bins_along_track', 'bins_across_track', 'number_layers',))
    sgm_dcol_h2o.units = 'molec./cm2'
    sgm_dcol_h2o.long_name = 'H2O layer column density'
    sgm_dcol_h2o.valid_min = 0.
    sgm_dcol_h2o.valid_max = 1.E+23
    sgm_dcol_h2o.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            sgm_dcol_h2o[ialt, iact, :] = atm[ialt][iact].H2O[:]

    # total coljmn of CO2, CH4, and H2O

    sgm_col_co2 = output_atm.createVariable('col_co2', np.float64, ('bins_along_track', 'bins_across_track',))
    sgm_col_co2.units = 'molec./cm2'
    sgm_col_co2.long_name = 'CO2 total column density'
    sgm_col_co2.valid_min = 0.
    sgm_col_co2.valid_max = 1.E+25
    sgm_col_co2.FillValue = -32767

    sgm_col_ch4 = output_atm.createVariable('col_ch4', np.float64, ('bins_along_track', 'bins_across_track',))
    sgm_col_ch4.units = 'molec./cm2'
    sgm_col_ch4.long_name = 'CH4 total column density'
    sgm_col_ch4.valid_min = 0.
    sgm_col_ch4.valid_max = 1.E+25
    sgm_col_ch4.FillValue = -32767

    sgm_col_h2o = output_atm.createVariable('col_h2o', np.float64, ('bins_along_track', 'bins_across_track',))
    sgm_col_h2o.units = 'molec./cm2'
    sgm_col_h2o.long_name = 'H2O total column density'
    sgm_col_h2o.valid_min = 0.
    sgm_col_h2o.valid_max = 1.E+25
    sgm_col_h2o.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            XAIR = np.sum(atm[ialt, iact].air[:])
            sgm_col_co2[ialt, iact] = np.sum(atm[ialt, iact].CO2[:])/XAIR*1.e6  # [ppmv]
            sgm_col_ch4[ialt, iact] = np.sum(atm[ialt, iact].CH4[:])/XAIR*1.e9  # [ppbv]
            sgm_col_h2o[ialt, iact] = np.sum(atm[ialt, iact].H2O[:])/XAIR*1.e6  # [ppmv]

    output_atm.close()

    return

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
    from end_to_end.lib import libSGM
    from end_to_end.lib import libATM
    from end_to_end.lib import libRT
    from end_to_end.lib import libSURF

    # first get the geometry data

    gm_data = get_gm_data(config['gm_input'])

    nact = gm_data['sza'][0].size
    nalt = len(gm_data['sza'])

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    albedo = np.zeros([nalt, nact])

    if(config['profile'] == 'individual_spectra'):
        albedo[0, :] = config['albedo'][:]

    if(config['profile'] == 'single_swath'):
        
        for iscen in range(config['numb_atm_scenes']+1):
            outofrange = (config['scene_trans_index'][iscen] > 100) & \
                (config['scene_trans_index'][iscen] < 0) 
            if(outofrange):
                sys.exit('config parameter scene_trans_index out of range')

        for iscen in range(config['numb_atm_scenes']):
             ind_start = config['scene_trans_index'][iscen]
             ind_end   = config['scene_trans_index'][iscen+1]            
             albedo[0, ind_start:ind_end] = config['albedo'][iscen]
        
    if((config['profile'] == 'S2_microHH')):

        # get collocated S2 data

        file_exists = os.path.isfile(config['S2_dump'])
        if(file_exists and (not(config['s2_forced']))):
            albedo = np.load(config['S2_dump'])
        else:
            albedo = libSGM.get_sentinel2_albedo(gm_data, config)
            np.save(config['S2_dump'], albedo)    # .npy extension is added if not given

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    if((config['profile'] == 'individual_spectra') or(config['profile'] == 'single_swath')):

        nlay = config['atmosphere']['nlay']  # number of layers
        dzlay = config['atmosphere']['dzlay']

        # we assume the same standard atmosphere for all pixels of the granule

        atm_std = libATM.get_AFGL_atm_homogenous_distribtution(config['afgl_input'], nlay, dzlay)

        atm = np.ndarray((nalt, nact), np.object_)
        for ialt in range(nalt):
            for iact in range(nact):
                atm[ialt, iact] = deepcopy(atm_std)

    if(config['profile'] == 'S2_microHH'):

        nlay = config['atmosphere']['nlay']  # number of layers
        dzlay = config['atmosphere']['dzlay']

        atm_std = libATM.get_AFGL_atm_homogenous_distribtution(config['afgl_input'], nlay, dzlay)

        if(config['only_afgl']):
            atm = atm_std
        else:         
            # get collocated mciroHH data

            if ((not os.path.exists(config['microHH_dump'])) or config['microHH_forced']):
                microHH = libATM.get_microHH_atm(gm_data['lat'], gm_data['lon'], config['microHH_data_path'],
                                                 config['microHH'],
                                                 config['kernel_parameter'])
                # Dump microHH dictionary into temporary pkl file
                pickle.dump(microHH, open(config['microHH_dump'], 'wb'))

            else:
                
                # Read microHH from pickle file
                microHH = pickle.load(open(config['microHH_dump'], 'rb'))

            atm = libATM.combine_microHH_standard_atm(microHH, atm_std)

    # =============================================================================
    #  Radiative transfer simulations
    # =============================================================================
    print('radiative transfer simuation...')

    # define line-by-line wavelength grid
    rad_output = {}
    rad_output['wavelength lbl'] = np.arange(config['spec_settings']['wave_start'], config['spec_settings']['wave_end'],
                                             config['spec_settings']['dwave'])  # nm

    nwav = len(rad_output['wavelength lbl'])
    # generate optics object for one representative model atmosphere of the domain

    nalt_ref = np.int0(nalt/2 - 0.5)
    nact_ref = np.int0(nact/2 - 0.5)

    optics = libRT.optic_abs_prop(rad_output['wavelength lbl'], atm[nalt_ref, nact_ref].zlay)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(rad_output['wavelength lbl'])

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
    rad_output['solar irradiance'] = np.interp(rad_output['wavelength lbl'], sun['wl'], sun['phsm2nm'])

    # Calculate surface data
    surface = libSURF.surface_prop(rad_output['wavelength lbl'])

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
    # sgm output to radiometric and geophysical output file
    # =============================================================================

    sgm_output(config['rad_output'], config['geo_output'], rad_output, atm, albedo)

    print('=>sgm calcultion finished successfully')
    return


if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    scene_generation_module(config)
