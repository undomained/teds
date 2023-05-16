# =============================================================================
# scene generation module for different E2E simulator profiles
# =============================================================================

import numpy as np
import sys
import os
import yaml
import pickle
import matplotlib.pyplot as plt
import netCDF4 as nc
from tqdm import tqdm
from copy import deepcopy

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


def sgm_output(path, filename_rad, filename_atm, rad_output, atm, albedo):

    nalt = len(rad_output['radiance'][:, 0, 0])
    nact = len(rad_output['radiance'][0, :, 0])
    nlbl = len(rad_output['radiance'][0, 0, :])

    file = path+filename_rad+'.nc'
    output_rad = nc.Dataset(file, mode='w')

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

    file = path+filename_atm+'.nc'
    output_atm = nc.Dataset(file, mode='w')

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

def scene_generation_module(paths, global_config, local_config):
    """

    Parameters
    ----------
    global_config : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from modules.lib import libSGM
    from modules.lib import libATM
    from modules.lib import libRT
    from modules.lib import libSURF
    from modules.lib import libNumTools

    run_id = '_'+global_config['run_id']

    # first get the geometry data

    gm_data_file = local_config['sgm_filename']['gm_input']+'_'+global_config['profile']+run_id

    gm_data = get_gm_data(paths.project+paths.data_interface+paths.interface_gm, gm_data_file)

    nact = gm_data['sza'][0].size
    nalt = len(gm_data['sza'])

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    albedo = np.empty([nalt, nact])

    if(global_config['profile'] == 'individual_spectra'):
        albedo[0, :] = global_config["individual_spectra"]["albedo"][:]

    if(global_config['profile'] == 'single_swath'):
        albedo[0, :] = global_config["single_swath"]["albedo"]

    if((global_config['profile'] == 'S2_hom_atm') or (global_config['profile'] == 'S2_microHH')):

        # get collocated S2 data
        albedo_tmp_file = 'sgm_albedo'+run_id+'.npy'

        file_exists = os.path.isfile(paths.project+paths.data_tmp+albedo_tmp_file)
        if(file_exists and (not(local_config['s2_forced']))):
            albedo = np.load(paths.project+paths.data_tmp+albedo_tmp_file)
        else:
            albedo = libSGM.get_sentinel2_albedo(gm_data, local_config)
            np.save(paths.project+paths.data_tmp+albedo_tmp_file, albedo)    # .npy extension is added if not given

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================

    if((global_config['profile'] == 'individual_spectra') or
       (global_config['profile'] == 'single_swath')):

        nlay = local_config['atmosphere']['nlay']  # number of layers
        dzlay = local_config['atmosphere']['dzlay']

        # we assume the same standard atmosphere for all pixels of the granule

        afgl_file = paths.project+paths.data_afgl+local_config['std_atm']
        atm_std = libATM.get_AFGL_atm_homogenous_distribtution(afgl_file, nlay, dzlay)

        atm = np.ndarray((nalt, nact), np.object_)
        for ialt in range(nalt):
            for iact in range(nact):
                atm[ialt, iact] = deepcopy(atm_std)

    if((global_config['profile'] == 'S2_microHH')):

        # get collocated mciroHH data
        microHH_tmp_file = 'microHH2_'+run_id+'.pkl'

        if ((not os.path.exists(paths.project + paths.data_tmp +microHH_tmp_file)) or local_config['microHH_forced']):
            data_path = paths.project+paths.data_microHH+ local_config['microHH']['sim']
            microHH = libATM.get_microHH_atm(gm_data['lat'], gm_data['lon'], data_path,
                                             local_config['microHH'],
                                             local_config['kernel_parameter'])
            # Dump microHH dictionary into temporary pkl file
            pickle.dump(microHH, open(paths.project+paths.data_tmp+microHH_tmp_file, 'wb'))

        else:

            # Read microHH from pickle file
            microHH = pickle.load(open(paths.project+paths.data_tmp+microHH_tmp_file, 'rb'))

        nlay = local_config['atmosphere']['nlay']  # number of layers
        dzlay = local_config['atmosphere']['dzlay']

        afgl_file = paths.project+paths.data_afgl+local_config['std_atm']
        atm_std = libATM.get_AFGL_atm_homogenous_distribtution(afgl_file, nlay, dzlay)

        atm = libATM.combine_microHH_standard_atm(microHH, atm_std)

        pickle.dump(atm, open(paths.project+paths.data_tmp+'microHH_plot.pkl', 'wb'))

    # =============================================================================
    #  Radiative transfer simulations
    # =============================================================================
    print('radiative transfer simuation...')
    # get cross sections
    xsec_file = paths.project + paths.data_tmp+'optics_prop'+run_id+'.pkl'

    # define line-by-line wavelength grid
    rad_output = {}
    rad_output['wavelength lbl'] = np.arange(local_config['spec_settings']['wave_start'], local_config['spec_settings']['wave_end'],
                                             local_config['spec_settings']['dwave'])  # nm

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
    if ((not os.path.exists(xsec_file)) or local_config['xsec_forced']):
        molec.get_data_HITRAN('./data/', iso_ids)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm[nalt_ref][nact_ref])
        # Dump optics.prop dictionary into temporary pkl file
        pickle.dump(optics.prop, open(xsec_file, 'wb'))
    else:
        # Read optics.prop dictionary from pickle file
        optics.prop = pickle.load(open(xsec_file, 'rb'))

    # solar irradiance spectrum
    sun = libRT.read_sun_spectrum_TSIS1HSRS(
        paths.project+ paths.data_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc')
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

    file_rad = local_config['sgm_filename']['filename_rad_output']  + '_' + global_config['profile']+run_id
    file_geo = local_config['sgm_filename']['filename_atm_output']  + '_' + global_config['profile']+run_id
    sgm_output(paths.project+paths.data_interface + paths.interface_sgm, file_rad, file_geo, rad_output, atm, albedo)

    print('=>sgm calcultion finished successfully for runid ' + global_config['run_id'])
    return
