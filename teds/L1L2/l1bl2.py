# ==============================================================================
#     level-1b to level-2 processor
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# ==============================================================================

import os
import sys
import numpy as np
import pickle as pkl
from tqdm import tqdm
from copy import deepcopy
import netCDF4 as nc
import yaml

from lib import libNumTools
from lib import libRT
from lib import libATM
from lib import libINV
from lib.libWrite import writevariablefromname
import matplotlib.pyplot as plt

class Emptyclass:
    """Empty class. Data container."""
    
    pass

def get_l1b(filename):
    # getting l1b data from file
    input = nc.Dataset(filename, mode='r')
    l1b_data = {}
    l1b_data['sza'] = deepcopy(input['GEOLOCATION_DATA']['sza'][:])
    l1b_data['saa'] = deepcopy(input['GEOLOCATION_DATA']['saa'][:])
    l1b_data['vza'] = deepcopy(input['GEOLOCATION_DATA']['vza'][:])
    l1b_data['vaa'] = deepcopy(input['GEOLOCATION_DATA']['vaa'][:])
    l1b_data['latitude'] = deepcopy(input['GEOLOCATION_DATA']['lat'][:])
    l1b_data['longitude'] = deepcopy(input['GEOLOCATION_DATA']['lon'][:])
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelength'][:])
    l1b_data['radiance'] = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
#    l1b_data['radmask'] = deepcopy(input['OBSERVATION_DATA']['radiance_mask'][:])
    l1b_data['noise'] = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])
    
    input.close()
    return (l1b_data)


def get_sgm_atm(filen_sgm_atm):
    data = nc.Dataset(filen_sgm_atm, mode='r')
    atm_sgm = {}
    surf_sgm = {}
    surf_sgm['albedo'] = deepcopy(data['albedo'][:])
    atm_sgm['zlay'] = deepcopy(data['zlay'][:])
    atm_sgm['zlev'] = deepcopy(data['zlev'][:])
    atm_sgm['dcol_co2'] = deepcopy(data['dcol_co2'][:])
    atm_sgm['dcol_ch4'] = deepcopy(data['dcol_ch4'][:])
    atm_sgm['dcol_h2o'] = deepcopy(data['dcol_h2o'][:])
    atm_sgm['col_co2'] = deepcopy(data['XCO2'][:])
    atm_sgm['col_ch4'] = deepcopy(data['XCH4'][:])
    atm_sgm['col_h2o'] = deepcopy(data['XH2O'][:])
    data.close()
    return (surf_sgm, atm_sgm)


def write_gasdata(gas, output_l2, _dims, _dims3d, nalt_l2, nact, nlay, l2product, retrieval_init):
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names
    xgas = "X" + gas
    xgas_precision = xgas + " precision"
    xgas_ker = xgas + ' col avg kernel'
    l2_X = np.zeros((nalt_l2, nact))
    l2_X_prec  = np.zeros((nalt_l2, nact))
    l2_X_avgk  = np.zeros((nalt_l2, nact, nlay))
    l2_X_prof  = np.zeros((nalt_l2, nact, nlay))

    for ialt in range(nalt_l2):
        for iact in range(nact):
            tmp = (l2product[ialt, iact][xgas]*scale[gas])
            l2_X[ialt, iact] = tmp
            l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
            l2_X_avgk[ialt, iact, :] = l2product[ialt, iact][xgas_ker][:]
            l2_X_prof[ialt, iact, :] = retrieval_init['trace gases'][gas]['ref_profile']
            
    # write data
    prec_varname   = 'precision'+xgas
    avgker_varname = 'avgkernel'+xgas
    aprof_varname  = 'apriori_profile'+gas

    _ = writevariablefromname(output_l2, xgas,           _dims,   l2_X)
    _ = writevariablefromname(output_l2, prec_varname,   _dims,   l2_X_prec)
    _ = writevariablefromname(output_l2, avgker_varname, _dims3d, l2_X_avgk)
    _ = writevariablefromname(output_l2, aprof_varname,  _dims3d, l2_X_prof)

def write_proxygasdata(gas, output_l2, _dims, nalt_l2, nact, nlay, l2product):
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names
    xgas = "X" + gas + " proxy"
    xgas_precision = xgas + " precision"
    l2_X = np.zeros((nalt_l2, nact))
    l2_X_prec = np.zeros((nalt_l2, nact))
    l2_X_accu = np.zeros((nalt_l2, nact))
    l2_X_qav  = np.zeros((nalt_l2, nact), dtype = np.int16)
    
    for ialt in range(nalt_l2):
        for iact in range(nact):
            tmp = (l2product[ialt, iact][xgas]*scale[gas])
            l2_X[ialt, iact] = tmp
            l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
            l2_X_accu[ialt, iact] = 99999.
            l2_X_qav[ialt, iact]  = 100

    # define names in constants_outputvariables
    vargas = "proxyX"+gas
    varprecgas = "precision"+vargas
    varaccugas = "accuracy" +vargas
    varqavgas  = "qa_value" +vargas
    # write data
    _ = writevariablefromname(output_l2, vargas, _dims, l2_X)
    _ = writevariablefromname(output_l2, varprecgas, _dims, l2_X_prec)
    _ = writevariablefromname(output_l2, varaccugas, _dims, l2_X_accu)
    _ = writevariablefromname(output_l2, varqavgas,  _dims, l2_X_qav)

def level2_output(filename, l2product, retrieval_init, l1bproduct, settings):
    # output of level 2 data
    nalt_l2 = len(l2product)
    nalt = l1bproduct['latitude'][:, 0].size
    nact = len(l2product[0])
    nlay = len(l2product[0][0]['XCO2 col avg kernel'])

    if (settings['sw_ALT_select']):
        nalt_start = settings['first_ALT_index']
        nalt_end = settings['last_ALT_index']
    else:
        nalt_start = 0
        nalt_end = nalt
    nalt_l2 = nalt_end-nalt_start
    # filename
    print(filename)
    output_l2 = nc.Dataset(filename, mode='w')
    output_l2.title = 'Tango Carbon E2ES L2 product'
    output_l2.createDimension('number_layers', nlay)     # spectral axis
    output_l2.createDimension('bins_across_track', nact)     # across track axis
    output_l2.createDimension('bins_along_track', nalt_l2)     # along track axis
    output_l2.createDimension('bins_albedo', 1)     # spectral bins albedo
    grp_ns    = output_l2.createGroup('non_scattering_retrieval')
    grp_prior = output_l2.createGroup('prior')
    grp_diag  = output_l2.createGroup('diagnostics')

    # layer height

    _ = writevariablefromname(output_l2, 'central_layer_height', ('number_layers',),  retrieval_init['zlay'])

    # dimensions
    _dimalt = ('bins_along_track')
    _dims   = ('bins_along_track', 'bins_across_track')
    _dims3d = ('bins_along_track', 'bins_across_track', 'number_layers',)
    # variables
    l2_conv = np.zeros((nalt_l2, nact), dtype = np.int32)
    l2_numb_iter = np.zeros((nalt_l2, nact), dtype = np.int32)
    l2_proc_flag = np.zeros((nalt_l2, nact), dtype = np.int16) + 100
    l2_chi2 = np.zeros((nalt_l2, nact))
    l2_alb = np.zeros((nalt_l2, nact))
    l2_spec_shift  = np.zeros((nalt_l2, nact))
    l2_spec_squeeze = np.zeros((nalt_l2, nact))
    l2_aqutime = np.zeros(nalt_l2) + 99999.   #aquisition time (dummy)
    
    for ialt in range(nalt_l2):
        for iact in range(nact):
            l2_conv[ialt, iact]         = int(l2product[ialt, iact]['convergence'])
            l2_numb_iter[ialt, iact]    = l2product[ialt, iact]['number_iter']
            l2_chi2[ialt, iact]         = l2product[ialt, iact]['chi2']
            l2_alb[ialt, iact]          = l2product[ialt, iact]['alb0']
            l2_spec_shift[ialt, iact]   = l2product[ialt, iact]['spec_shift']
            l2_spec_squeeze[ialt, iact] = l2product[ialt, iact]['spec_squeeze']

    # aquisition time
    _ = writevariablefromname(output_l2, 'aqui_time',       _dimalt, l2_aqutime)
    # processing flag
    _ = writevariablefromname(grp_diag,  'process_flag',    _dims, l2_proc_flag)
    # convergence
    _ = writevariablefromname(grp_diag,  'convergence',     _dims, l2_conv)
    # latitude
    _ = writevariablefromname(output_l2, 'latitude',        _dims, l1bproduct['latitude'][nalt_start:nalt_end, :nact])
    # longitude
    _ = writevariablefromname(output_l2, 'longitude',       _dims, l1bproduct['longitude'][nalt_start:nalt_end, :nact])
    # maxiterations
    _ = writevariablefromname(grp_diag,  'maxiterations',   _dims, l2_numb_iter)
    # chi2
    _ = writevariablefromname(grp_diag,  'spectralchi2',    _dims, l2_chi2)
    # albedo
    _ = writevariablefromname(output_l2, 'albedo',          _dims, l2_alb)
    # surface pressure
    _ = writevariablefromname(grp_prior, 'surface_pressure',_dims, retrieval_init['surface pressure'])
    # surface elevation
    _ = writevariablefromname(grp_prior, 'surface_elevation',_dims, retrieval_init['surface elevation'])
    # spectral shift
    _ = writevariablefromname(output_l2, 'spectralshift',   _dims, l2_spec_shift)
    # chi2
    _ = writevariablefromname(output_l2, 'spectralsqueeze', _dims, l2_spec_squeeze)
    
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names

    gases = ['CO2', 'CH4', 'H2O']
    for gas in gases:
        xgas = "X" + gas
        xgas_precision = xgas + " precision"
        xgas_ker = xgas + ' col avg kernel'

        l2_X = np.zeros((nalt_l2, nact))
        l2_X_prec  = np.zeros((nalt_l2, nact))
        l2_X_avgk  = np.zeros((nalt_l2, nact, nlay))
        l2_X_prof  = np.zeros((nalt_l2, nact, nlay))

        for ialt in range(nalt_l2):
            for iact in range(nact):
                tmp = (l2product[ialt, iact][xgas]*scale[gas])
                l2_X[ialt, iact] = tmp
                l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
                l2_X_avgk[ialt, iact, :] = l2product[ialt, iact][xgas_ker][:]
                l2_X_prof[ialt, iact, :] = retrieval_init['trace gases'][gas]['ref_profile']
            
            # write data
            prec_varname   = 'precision'+xgas
            avgker_varname = 'avgkernel'+xgas
            aprof_varname  = 'apriori_profile'+gas

        _ = writevariablefromname(grp_ns, xgas,           _dims,   l2_X)
        _ = writevariablefromname(grp_ns, prec_varname,   _dims,   l2_X_prec)
        _ = writevariablefromname(output_l2, avgker_varname, _dims3d, l2_X_avgk)
        _ = writevariablefromname(grp_prior, aprof_varname,  _dims3d, l2_X_prof)
    
    #next the main product, proxy
    gases = ['CO2', 'CH4']
    for gas in gases:

        xgas = "X" + gas + " proxy"
        xgas_precision = xgas + " precision"
        l2_X = np.zeros((nalt_l2, nact))
        l2_X_prec = np.zeros((nalt_l2, nact))
        l2_X_accu = np.zeros((nalt_l2, nact))
        l2_X_qav  = np.zeros((nalt_l2, nact), dtype = np.int16)
    
        for ialt in range(nalt_l2):
            for iact in range(nact):
                l2_X[ialt, iact] = (l2product[ialt, iact][xgas]*scale[gas])
                l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
                l2_X_accu[ialt, iact] = 99999.
                l2_X_qav[ialt, iact]  = 100

        # define names in constants_outputvariables
        vargas = "proxyX"+gas
        varprecgas = "precision"+vargas
        varaccugas = "accuracy" +vargas
        varqavgas  = "qa_value" +vargas
        # write data
        _ = writevariablefromname(output_l2, vargas, _dims, l2_X)
        _ = writevariablefromname(output_l2, varprecgas, _dims, l2_X_prec)
        _ = writevariablefromname(output_l2, varaccugas, _dims, l2_X_accu)
        _ = writevariablefromname(output_l2, varqavgas,  _dims, l2_X_qav)

    output_l2.close()


def level2_diags_output(filename, l2product, measurement):
    # this function writes some ectra diagnostics to a file that are not included in the level 2 data product.
    nalt = len(l2product)
    nact = len(l2product[0])
    nwave = len(measurement[0, 0]['wavelength'])
    # file name
    print('diag_output')
    output_l2diag = nc.Dataset(filename, mode='w')
    output_l2diag.title = 'Tango Carbon E2ES L2 diagnostics'
    output_l2diag.createDimension('bins_spectral', nwave)     # spectral axis
    output_l2diag.createDimension('bins_along_track', nalt)     # across track axis
    output_l2diag.createDimension('bins_across_track', nact)     # across track axis

    # create 2d arrays
    l2diag_wave = np.zeros((nalt, nact, nwave))
    l2diag_measurement = np.zeros((nalt, nact, nwave))
    l2diag_sun_irradiance = np.zeros((nalt, nact, nwave))
    l2diag_gainCO2 = np.zeros((nalt, nact, nwave))
    l2diag_gainH2O = np.zeros((nalt, nact, nwave))
    l2diag_gainCH4 = np.zeros((nalt, nact, nwave))
    for ialt in range(nalt):
        for iact in range(nact):
            l2diag_wave[ialt, iact, :] = measurement[ialt, iact]['wavelength'][:]
            l2diag_measurement[ialt, iact, :] = measurement[ialt, iact]['ymeas'][:]
            l2diag_sun_irradiance[ialt, iact, :] = measurement[ialt, iact]['sun'][:]
            l2diag_gainCO2[ialt, iact, :] = l2product[ialt, iact]['gain_XCO2'][:]
            l2diag_gainH2O[ialt, iact, :] = l2product[ialt, iact]['gain_XH2O'][:]
            l2diag_gainCH4[ialt, iact, :] = l2product[ialt, iact]['gain_XCH4'][:]

    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    # wavelength of measurement
    _ = writevariablefromname(output_l2diag, "l2wavelength", _dims, l2diag_wave)
    # spectral radiance measurements
    _ = writevariablefromname(output_l2diag, "l2measurement", _dims, l2diag_measurement)
    # solar irradiance spectrum
    _ = writevariablefromname(output_l2diag, "l2solarirradiance", _dims, l2diag_sun_irradiance)
    # gain CO2
    _ = writevariablefromname(output_l2diag, "l2gainCO2", _dims, l2diag_gainCO2)
    # gain CH4
    _ = writevariablefromname(output_l2diag, "l2gainCH4", _dims, l2diag_gainCH4)
    # gain H2O
    _ = writevariablefromname(output_l2diag, "l2gainH2O", _dims, l2diag_gainH2O)
    output_l2diag.close()

def level1b_to_level2_processor(config, sw_diag_output = False):
    # get the l1b files

    print('level 1B to 2 proessor ...')
    l1b = get_l1b(config['l1b_input'])

    # get sgm geo data
    surf_sgm, atm_sgm = get_sgm_atm(config['sgm_input'])

    # get pixel mask
    # if (config['retrieval_init']['sw_pixel_mask']):
#        print('take pixel mask')get_AFGL_atm_homogenous_distribution
#        mask = np.load(config['pixel_mask_input'])
# 
    # else:
        # nact = l1b['radiance'][0, :, 0].size
        # nwave = l1b['radiance'][0, 0, :].size
        # mask = np.full((nact, nwave), True)

    #mask on radiance nan values
    mask = np.logical_not(np.isnan(l1b['radiance']))

    # fig = plt.figure(figsize=(10, 8), dpi=100)
    # radiance = l1b['radiance'][0, :,:]
    # masked_data = np.ma.masked_where( np.invert(mask).data, radiance)
    # plt.imshow(masked_data, vmin = 1.E16,vmax = 4.E16, cmap = 'viridis', aspect = 4)

    # Internal lbl spectral grid
    wave_start = config['spec_settings']['wavestart']
    wave_end = config['spec_settings']['waveend']
    wave_extend = config['spec_settings']['wave_extend']
    dwave_lbl = config['spec_settings']['dwave']
    wave_lbl = np.arange(wave_start-wave_extend, wave_end+wave_extend, dwave_lbl)  # nm

    # Vertical layering
    nlay = config['atmosphere']['nlay']
    dzlay = config['atmosphere']['dzlay']
    psurf = config['atmosphere']['psurf']
    nlev = nlay + 1  # number of levels

    zlay = (np.arange(nlay-1, -1, -1)+0.5)*dzlay  # altitude of layer midpoint
    zlev = np.arange(nlev-1, -1, -1)*dzlay  # altitude of layer interfaces = levels

    # model atmosphere

    atm = libATM.atmosphere_data(zlay, zlev, psurf)
    atm.get_data_AFGL(config['afgl_input'])

    # XCO2 = np.sum(atm.CO2)/np.sum(atm.air)
    # XCH4 = np.sum(atm.CH4)/np.sum(atm.air)
    # XH2O = np.sum(atm.H2O)/np.sum(atm.air)
    # print(XCO2*1.E6, XCH4*1.E9, XH2O*1E6)

    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['sun_reference'])
    sun_lbl = libRT.interpolate_sun(sun, wave_lbl)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(wave_lbl)
    molec.get_data_HITRAN(config['hapi_path'], iso_ids)

    # Calculate optical properties
    # If pickle file exists read from file
    if ((not os.path.exists(config['xsec_dump'])) or config['xsec_forced']):
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm)
        # Dump optics.prop dictionary into temporary pkl file
        pkl.dump(optics.prop, open(config['xsec_dump'], 'wb'))
    else:
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Read optics.prop dictionary from pickle file
        optics.prop = pkl.load(open(config['xsec_dump'], 'rb'))

    optics.set_opt_depth_species(atm, ['molec_01', 'molec_32', 'molec_07'])

    # Initialization of the least squares fit
    retrieval_init = {}
    retrieval_init['chi2 limit'] = config['retrieval_init']['chi2_lim']
    retrieval_init['maximum iteration'] = config['retrieval_init']['max_iter']
    retrieval_init['zlay'] = zlay
    retrieval_init['trace gases'] = {'CO2': {'init': 300,      'scaling': 1E-6},
                                     'CH4': {'init': 1700,     'scaling': 1E-9},
                                     'H2O': {'init': 1000,     'scaling': 1E-6}}

#    retrieval_init['surface'] = {'alb0': 0.17}
    retrieval_init['wavelength lbl'] = wave_lbl
    retrieval_init['solar irradiance'] = sun_lbl

    nact = l1b['latitude'][0, :].size
    nalt = l1b['latitude'][:, 0].size

    l2product = {}
    xch4_model = 1.8E-6
    xco2_model = 405.E-6

    XCO2 = np.zeros(nact)
    XCO2_true_smoothed = np.zeros(nact)
    XCO2_true = np.zeros(nact)
    XCO2_prec = np.zeros(nact)

    if (config['retrieval_init']['sw_ALT_select']):
        nalt_start = config['retrieval_init']['first_ALT_index']
        nalt_end = config['retrieval_init']['last_ALT_index']
    else:
        nalt_start = 0
        nalt_end = nalt
    nalt_l2 = nalt_end-nalt_start
    
    l2product = np.ndarray((nalt_l2, nact), np.object_)
    measurement = np.ndarray((nalt_l2, nact), np.object_)

    # plt.plot(l1b['wavelength'][0,mask[0,:]], l1b['radiance'][0,0,mask[0,:]])
    # plt.plot(l1b['wavelength'][10,mask[10,:]], l1b['radiance'][0,10,mask[10,:]])
    # sys.exit()
    # We introduce two types of ialt indices, ialt points to l1b data structure and ilat_l2 to l2 data structure.
    # Later is different to l1b structure because of option for image selection (sw_ALT_select).

    runtime_cum = {}
    runtime_cum['opt'] = 0.
    runtime_cum['rtm'] = 0.
    runtime_cum['conv'] = 0.
    runtime_cum['kern'] = 0.

    for ialt_l2, ialt in enumerate(tqdm(range(nalt_start, nalt_end))):
        for iact in range(nact):
            # number of 'good' pixels
            numb_spec_points = np.sum(mask[ialt_l2,iact,:])/float(mask[ialt_l2,iact,:].size)
            if(numb_spec_points >=0.9):
                # initialization of pixel retrieval
                xco2_ref = 405.  # ppm
                xco2 = np.sum(atm.CO2) / np.sum(atm.air) * 1.E6
    
                if (config['retrieval_init']['sw_prof_init'] == 'afgl'):
                    retrieval_init['trace gases']['CO2']['ref_profile'] = atm.CO2*xco2_ref/xco2
                    retrieval_init['trace gases']['CH4']['ref_profile'] = atm.CH4
                    retrieval_init['trace gases']['H2O']['ref_profile'] = atm.H2O
                elif (config['retrieval_init']['sw_prof_init'] == 'sgm'):
                    retrieval_init['trace gases']['CO2']['ref_profile'] = atm_sgm['dcol_co2'][ialt, iact, :]*xco2_ref/xco2
                    retrieval_init['trace gases']['CH4']['ref_profile'] = atm_sgm['dcol_ch4'][ialt, iact, :]
                    retrieval_init['trace gases']['H2O']['ref_profile'] = atm_sgm['dcol_h2o'][ialt, iact, :]
                else:
                    sys.exit('sw_prof_init of l1bl2 configuration not set correctly!')
                    
                retrieval_init['surface pressure']  = np.zeros([nalt_l2,nact])+1013.  #these are dummy values for the time being
                retrieval_init['surface elevation'] = np.zeros([nalt_l2,nact])

                wavelength = l1b['wavelength'][iact, mask[ialt,iact, :]].data
                istart = np.argmin(np.abs(wavelength - wave_start))
                iend = np.argmin(np.abs(wavelength - wave_end))
                wave_meas = wavelength[istart:iend+1]  # nm
    
                # define isrf function
    
                isrf_convolution = libNumTools.get_isrf(wave_meas, wave_lbl, config['isrf_settings'])
    
                atm_ret = deepcopy(atm)  # to initialize each retrieval with the same atmosphere
    
                sun = isrf_convolution(sun_lbl)
    
                # Observation geometry
                nwave = wave_meas.size
                measurement[ialt_l2, iact] = {}
                measurement[ialt_l2, iact]['wavelength'] = wave_meas[:]
                measurement[ialt_l2, iact]['mu0'] = np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
                measurement[ialt_l2, iact]['muv'] = np.cos(np.deg2rad(l1b['vza'][ialt, iact]))
                measurement[ialt_l2, iact]['ymeas'] = l1b['radiance'][ialt, iact, mask[ialt,iact, :]][istart:iend+1].data
                measurement[ialt_l2, iact]['Smeas'] = np.eye(
                    nwave)*(l1b['noise'][ialt, iact, mask[ialt,iact, :]][istart:iend+1].data)**2
                measurement[ialt_l2, iact]['sun'] = sun
    
                # derive first guess albedo from the maximum reflectance
                ymeas_max = np.max(measurement[ialt_l2, iact]['ymeas'])
    
                idx = np.where(measurement[ialt_l2, iact]['ymeas'] == ymeas_max)[0][0]
                alb_first_guess = measurement[ialt_l2, iact]['ymeas'][idx] / \
                    sun[idx]*np.pi/np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
                retrieval_init['surface'] = {'alb0': alb_first_guess, 'alb1': 0.0}
    
                # Non-scattering least squares fit
                l2product[ialt_l2, iact], runtime = libINV.Gauss_Newton_iteration(
                    retrieval_init, atm_ret, optics, measurement[ialt_l2, iact],
                    config['retrieval_init']['max_iter'], config['retrieval_init']['chi2_lim'],
                    isrf_convolution)
    
                for key in runtime.keys():
                    runtime_cum[key] = runtime_cum[key] + runtime[key]
    
                if (not l2product[ialt_l2, iact]['convergence']):
                    print('pixel did not converge (ialt,iact) = ', ialt_l2, iact)
    
                # Define proxy product
                l2product[ialt_l2, iact]['XCO2 proxy'] = l2product[ialt_l2,
                                                                   iact]['XCO2']/l2product[ialt_l2, iact]['XCH4']*xch4_model
                l2product[ialt_l2, iact]['XCH4 proxy'] = l2product[ialt_l2,
                                                                   iact]['XCH4']/l2product[ialt_l2, iact]['XCO2']*xco2_model
                rel_error = np.sqrt((l2product[ialt_l2, iact]['XCO2 precision']/l2product[ialt_l2, iact]['XCO2'])**2 +
                                    (l2product[ialt_l2, iact]['XCH4 precision']/l2product[ialt_l2, iact]['XCH4'])**2)
                l2product[ialt_l2, iact]['XCO2 proxy precision'] = rel_error * l2product[ialt_l2, iact]['XCO2']
                l2product[ialt_l2, iact]['XCH4 proxy precision'] = rel_error * l2product[ialt_l2, iact]['XCH4']

            else:
                l2product[ialt_l2, iact] = level2_nan(retrieval_init, nlay)
                l2product[ialt_l2, iact]['XCO2 proxy'] = float("nan")
                l2product[ialt_l2, iact]['XCH4 proxy'] = float("nan")
                l2product[ialt_l2, iact]['XCO2 proxy precision'] = float("nan")
                l2product[ialt_l2, iact]['XCH4 proxy precision'] = float("nan")

            l2product[ialt_l2, iact]['spec_shift']   = 0.  #dummies for the time beeing
            l2product[ialt_l2, iact]['spec_squeeze'] = 0.

            # XCO2[iact] = l2product[ialt_l2, iact]['XCO2 proxy']*1.E6
            # XCO2_prec[iact] = l2product[ialt_l2, iact]['XCO2 proxy precision']*1.E6
            # XCO2_true_smoothed[iact] = np.dot(l2product[ialt_l2, iact]['XCO2 col avg kernel'],
            #                                   atm_sgm['dcol_co2'][ialt, iact, :])/np.sum(atm.air)*1.E6
            # XCO2_true[iact] = np.sum(atm_sgm['dcol_co2'][ialt, iact, :])/np.sum(atm.air)*1.E6

    # output to netcdf file
    level2_output(config['l2_output'], l2product, retrieval_init, l1b, config['retrieval_init'])

    if(sw_diag_output):
        level2_diags_output(config['l2_diag'], l2product, measurement)

    print('=> l1bl2 finished successfully')

#    print('cumulative run time: ',runtime_cum)

    return l2product

def level2_nan(retrieval_init, nlay):

    output = {}
    output['chi2'] = float("nan")
    output['convergence'] = False
    output['number_iter'] = 0.

    for l, key in enumerate(retrieval_init['trace gases'].keys()):
        output['X'+key] = float("nan")
        output['X'+key+' precision'] = float("nan")
        output['gain_X'+key] = float("nan")

    m = 0
    while 'alb%d' % (m) in retrieval_init['surface'].keys():
        output['alb%d' % (m)] = float("nan")
        m = m+1
    output['XCO2 col avg kernel'] = np.zeros(nlay)
    output['XCH4 col avg kernel'] = np.zeros(nlay)
    output['XH2O col avg kernel'] = np.zeros(nlay)

    return(output)
    
if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    level1b_to_level2_processor(config)
