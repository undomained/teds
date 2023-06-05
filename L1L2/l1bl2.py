import os
import sys
import numpy as np
import pickle as pkl
from tqdm import tqdm
from copy import deepcopy
import netCDF4 as nc

from end_to_end.lib import libNumTools
from end_to_end.lib import libRT
from end_to_end.lib import libATM
from end_to_end.lib import libINV
from end_to_end.lib import libNumTools
import matplotlib.pyplot as plt

def get_l1b(path, filename):

    #getting l1b data from file
    
    file = path+filename+'.nc'

    input = nc.Dataset(file, mode='r')

    l1b_data = {}
    l1b_data['sza'] = deepcopy(input['GEOLOCATION_DATA']['sza'][:])
    l1b_data['saa'] = deepcopy(input['GEOLOCATION_DATA']['saa'][:])
    l1b_data['vza'] = deepcopy(input['GEOLOCATION_DATA']['vza'][:])
    l1b_data['vaa'] = deepcopy(input['GEOLOCATION_DATA']['vaa'][:])
    l1b_data['latitude'] = deepcopy(input['GEOLOCATION_DATA']['lat'][:])
    l1b_data['longitude'] = deepcopy(input['GEOLOCATION_DATA']['lon'][:])
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelengths'][:])
    l1b_data['radiance'] = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    l1b_data['noise'] = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])

    input.close()

    return(l1b_data)

def level2_output(path, filename, l2product, retrieval_init, l1bproduct):

    #output of level 2 data
    
    nalt = len(l2product)
    nact = len(l2product[0])
    nlay = len(l2product[0][0]['XCO2 col avg kernel'])

    file = path+filename+'.nc'
    output_l2 = nc.Dataset(file, mode='w')

    output_l2.title = 'Tango Carbon E2ES L2 product'
    output_l2.createDimension('number_layers', nlay)     # spectral axis
    output_l2.createDimension('bins_across_track', nact)     # across track axis
    output_l2.createDimension('bins_along_track', nalt)     # along track axis
                      
    l2_conv = output_l2.createVariable('convergence', np.int32, ('bins_along_track', 'bins_across_track',))
    l2_conv.units = '1'
    l2_conv.long_name = 'convergence'
    l2_conv.valid_min = 0
    l2_conv.valid_max = 1
    l2_conv.FillValue = -32767
    for ialt in range(nalt):
        for iact in range(nact):
            l2_conv[ialt, iact] = int(l2product[ialt, iact]['convergence'])

    l2_zlay = output_l2.createVariable('zlay', np.float64, ('number_layers',))
    l2_zlay.units = 'm'
    l2_zlay.long_name = 'layer height'
    l2_zlay.valid_min = 0.
    l2_zlay.valid_max = 1.E6
    l2_zlay.FillValue = -32767
    l2_zlay[:] = retrieval_init['zlay'][:]

    l2_lat = output_l2.createVariable('latitude', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_lat.units = 'degree'
    l2_lat.long_name = 'latitude'
    l2_lat.valid_min = -90.
    l2_lat.valid_max = +90.
    l2_lat.FillValue = -32767
    l2_lat[:, :] = l1bproduct['latitude'][:nalt, :nact]

    l2_lon = output_l2.createVariable('longitude', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_lon.units = 'degree'
    l2_lon.long_name = 'longitude'
    l2_lon.valid_min = -90.
    l2_lon.valid_max = +90.
    l2_lon.FillValue = -32767
    l2_lon[:, :] = l1bproduct['longitude'][:nalt, :nact]

    l2_numb_iter = output_l2.createVariable('max iter', np.int32, ('bins_along_track', 'bins_across_track',))
    l2_numb_iter.units = '1'
    l2_numb_iter.long_name = 'number of iterations'
    l2_numb_iter.valid_min = 0
    l2_numb_iter.valid_max = retrieval_init['maximum iteration']
    l2_numb_iter.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2_numb_iter[ialt, iact] = l2product[ialt, iact]['number_iter']

    l2_chi2 = output_l2.createVariable('chi2', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_chi2.units = '1'
    l2_chi2.long_name = 'spectral chi square of the fit'
    l2_chi2.valid_min = 0.
    l2_chi2.valid_max = 100.
    l2_chi2.FillValue = -32767
    for ialt in range(nalt):
        for iact in range(nact):
            l2_chi2[ialt, iact] = l2product[ialt, iact]['chi2']

    l2_alb = output_l2.createVariable('albedo', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_alb.units = '1'
    l2_alb.long_name = 'Lambertian surface albedo'
    l2_alb.valid_min = 0.
    l2_alb.valid_max = 1.
    l2_alb.FillValue = -32767
    for ialt in range(nalt):
        for iact in range(nact):
            l2_alb[ialt, iact] = l2product[ialt, iact]['alb0']

    units = {'CO2': 'ppm', 'CH4': 'ppb', 'H2O': 'ppm'}
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    maxval = {'CO2': 2.E3, 'CH4': 8.E3, 'H2O': 5.E6}
    # non-scattering data product

    l2_XCO2 = output_l2.createVariable('XCO2', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCO2.units = units['CO2']
    l2_XCO2.long_name = 'CO2 dry air column mixing ratio XCO2'
    l2_XCO2.valid_min = 0.
    l2_XCO2.valid_max = 2.E20  # maxval['CO2']
    l2_XCO2.FillValue = -32767

    l2_XCO2_prec = output_l2.createVariable('precision XCO2', np.float64,
                                            ('bins_along_track', 'bins_across_track',))
    l2_XCO2_prec.units = units['CO2']
    l2_XCO2_prec.long_name = 'CO2 precision of dry air column mixing ratio XCO2'
    l2_XCO2_prec.valid_min = 0.
    l2_XCO2_prec.valid_max = 0.1*maxval['CO2']
    l2_XCO2_prec.FillValue = -32767

    l2_XCO2_avgk = output_l2.createVariable('col avg kernel XCO2', np.float64,
                                            ('bins_along_track', 'bins_across_track', 'number_layers',))
    l2_XCO2_avgk.units = '1'
    l2_XCO2_avgk.long_name = 'CO2 column averaging kernel of XCO2'
    l2_XCO2_avgk.valid_min = 0.
    l2_XCO2_avgk.valid_max = 10.
    l2_XCO2_avgk.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            tmp = (l2product[ialt, iact]['XCO2']*scale['CO2'])
            l2_XCO2[ialt, iact] = tmp
            l2_XCO2_prec[ialt, iact] = l2product[ialt, iact]['XCO2 precision']*scale['CO2']
            l2_XCO2_avgk[ialt, iact, :] = l2product[ialt, iact]['XCO2 col avg kernel'][:]

    l2_XCH4 = output_l2.createVariable('XCH4', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCH4.units = units['CH4']
    l2_XCH4.long_name = 'CH4 dry air column mixing ratio XCH4'
    l2_XCH4.valid_min = 0.
    l2_XCH4.valid_max = maxval['CH4']
    l2_XCH4.FillValue = -32767

    l2_XCH4_prec = output_l2.createVariable('precision XCH4', np.float64,
                                            ('bins_along_track', 'bins_across_track',))
    l2_XCH4_prec.units = units['CH4']
    l2_XCH4_prec.long_name = 'CH4 precision of dry air column mixing ratio XCH4'
    l2_XCH4_prec.valid_min = 0.
    l2_XCH4_prec.valid_max = 0.1*maxval['CH4']
    l2_XCH4_prec.FillValue = -32767

    l2_XCH4_avgk = output_l2.createVariable('col avg kernel XCH4', np.float64,
                                            ('bins_along_track', 'bins_across_track', 'number_layers',))
    l2_XCH4_avgk.units = '1'
    l2_XCH4_avgk.long_name = 'CH4 column averaging kernel of XCH4'
    l2_XCH4_avgk.valid_min = 0.
    l2_XCH4_avgk.valid_max = 10.
    l2_XCH4_avgk.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2_XCH4[ialt, iact] = (l2product[ialt, iact]['XCH4']*scale['CH4'])
            l2_XCH4_prec[ialt, iact] = l2product[ialt, iact]['XCH4 precision']*scale['CH4']
            l2_XCH4_avgk[ialt, iact, :] = l2product[ialt, iact]['XCH4 col avg kernel'][:]

    l2_XH2O = output_l2.createVariable('XH2O', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XH2O.units = units['H2O']
    l2_XH2O.long_name = 'H2O dry air column mixing ratio XH2O'
    l2_XH2O.valid_min = 0.
    l2_XH2O.valid_max = maxval['H2O']
    l2_XH2O.FillValue = -32767

    l2_XH2O_prec = output_l2.createVariable('precision XH2O', np.float64,
                                            ('bins_along_track', 'bins_across_track',))
    l2_XH2O_prec.units = units['H2O']
    l2_XH2O_prec.long_name = 'H2O precision of dry air column mixing ratio XH2O'
    l2_XH2O_prec.valid_min = 0.
    l2_XH2O_prec.valid_max = 0.1*maxval['H2O']
    l2_XH2O_prec.FillValue = -32767

    l2_XH2O_avgk = output_l2.createVariable('col avg kernel XH2O', np.float64,
                                            ('bins_along_track', 'bins_across_track', 'number_layers',))
    l2_XH2O_avgk.units = '1'
    l2_XH2O_avgk.long_name = 'H2O column averaging kernel of XH2O'
    l2_XH2O_avgk.valid_min = 0.
    l2_XH2O_avgk.valid_max = 10.
    l2_XH2O_avgk.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2_XH2O[ialt, iact] = l2product[ialt, iact]['XH2O']*scale['H2O']
            l2_XH2O_prec[ialt, iact] = l2product[ialt, iact]['XH2O precision']*scale['H2O']
            l2_XH2O_avgk[ialt, iact, :] = l2product[ialt, iact]['XH2O col avg kernel'][:]

    # proxy data product

    l2_XCO2_proxy = output_l2.createVariable('XCO2 proxy', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCO2_proxy.units = units['CO2']
    l2_XCO2_proxy.long_name = 'CO2 proxy dry air column mixing ratio XCO2'
    l2_XCO2_proxy.valid_min = 0.
    l2_XCO2_proxy.valid_max = maxval['CO2']
    l2_XCO2_proxy.FillValue = -32767

    l2_XCO2_proxy_prec = output_l2.createVariable('precision XCO2 proxy',
                                                  np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCO2_proxy_prec.units = units['CO2']
    l2_XCO2_proxy_prec.long_name = 'CO2 precision of proxy  dry air column mixing ratio XCO2'
    l2_XCO2_proxy_prec.valid_min = 0.
    l2_XCO2_proxy_prec.valid_max = 0.1*maxval['CO2']
    l2_XCO2_proxy_prec.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2_XCO2_proxy[ialt, iact] = l2product[ialt, iact]['XCO2 proxy']*scale['CO2']
            l2_XCO2_proxy_prec[ialt, iact] = l2product[ialt, iact]['XCO2 proxy precision']*scale['CO2']

    l2_XCH4_proxy = output_l2.createVariable('XCH4 proxy', np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCH4_proxy.units = units['CH4']
    l2_XCH4_proxy.long_name = 'CH4 proxy dry air column mixing ratio XCH4'
    l2_XCH4_proxy.valid_min = 0.
    l2_XCH4_proxy.valid_max = maxval['CH4']
    l2_XCH4_proxy.FillValue = -32767

    l2_XCH4_proxy_prec = output_l2.createVariable('precision XCH4 proxy',
                                                  np.float64, ('bins_along_track', 'bins_across_track',))
    l2_XCH4_proxy_prec.units = units['CH4']
    l2_XCH4_proxy_prec.long_name = 'CH4 precision of proxy  dry air column mixing ratio XCH4'
    l2_XCH4_proxy_prec.valid_min = 0.
    l2_XCH4_proxy_prec.valid_max = 0.1*maxval['CH4']
    l2_XCH4_proxy_prec.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2_XCH4_proxy[ialt, iact] = l2product[ialt, iact]['XCH4 proxy']*scale['CH4']
            l2_XCH4_proxy_prec[ialt, iact] = l2product[ialt, iact]['XCH4 proxy precision']*scale['CH4']

    output_l2.close()
    return


def level2_diags_output(path, filename, l2product, measurement):

    #this function writes some ectra diagnostics to a file that are not included in the level 2 data product.
    
    nalt = len(l2product)
    nact = len(l2product[0])
    nwave = len(measurement[0, 0]['wavelength'])

    file = path+filename+'.nc'
    output_l2diag = nc.Dataset(file, mode='w')

    output_l2diag.title = 'Tango Carbon E2ES L2 diagnostics'
    output_l2diag.createDimension('bins_spectral', nwave)     # spectral axis
    output_l2diag.createDimension('bins_along_track', nalt)     # across track axis
    output_l2diag.createDimension('bins_across_track', nact)     # across track axis

    l2diag_wave = output_l2diag.createVariable(
        'wavelength', np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral',))
    l2diag_wave.units = 'nm'
    l2diag_wave.long_name = 'wavelength of measurement'
    l2diag_wave.valid_min = 0.
    l2diag_wave.valid_max = 4000.
    l2diag_wave.FillValue = -32767
    
    l2diag_measurement = output_l2diag.createVariable('measurement',
                                                      np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral'))
    l2diag_measurement.units = 'photons/(nm m2 s sr)'
    l2diag_measurement.long_name = 'spetral radiance measurements'
    l2diag_measurement.valid_min = 0.
    l2diag_measurement.valid_max = 1.E26
    l2diag_measurement.FillValue = -32767

    l2diag_sun_irradiance = output_l2diag.createVariable('solar irradiance',
                                                      np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral'))
    l2diag_sun_irradiance.units = 'photons/(nm m2 s)'
    l2diag_sun_irradiance.long_name = 'solar irradiance spectrum'
    l2diag_sun_irradiance.valid_min = 0.
    l2diag_sun_irradiance.valid_max = 1.E26
    l2diag_sun_irradiance.FillValue = -32767

    l2diag_gainCO2 = output_l2diag.createVariable('gain CO2',
                                                  np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral'))
    l2diag_gainCO2.units = 'ppm/(photons/(nm m2 s sr))'
    l2diag_gainCO2.long_name = 'CO2 spectral gain vector'
    l2diag_gainCO2.valid_min = -1.E30
    l2diag_gainCO2.valid_max = 1.E30
    l2diag_gainCO2.FillValue = -32767

    l2diag_gainCH4 = output_l2diag.createVariable('gain CH4',
                                                  np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral'))
    l2diag_gainCH4.units = 'ppb/(photons/(nm m2 s sr))'
    l2diag_gainCH4.long_name = 'CH4 spectral gain vector'
    l2diag_gainCH4.valid_min = -1.E30
    l2diag_gainCH4.valid_max = +1.E30
    l2diag_gainCH4.FillValue = -32767

    l2diag_gainH2O = output_l2diag.createVariable('gain H2O',
                                                  np.float64, ('bins_along_track', 'bins_across_track', 'bins_spectral'))
    l2diag_gainH2O.units = 'ppm/(photons/(nm m2 s sr))'
    l2diag_gainH2O.long_name = 'CO2 spectral gain vector'
    l2diag_gainH2O.valid_min = -1.E30
    l2diag_gainH2O.valid_max = +1.E30
    l2diag_gainH2O.FillValue = -32767

    for ialt in range(nalt):
        for iact in range(nact):
            l2diag_wave[ialt, iact, :] = measurement[ialt, iact]['wavelength'][:]
            l2diag_measurement[ialt, iact, :] = measurement[ialt, iact]['ymeas'][:]
            l2diag_sun_irradiance[ialt, iact, :] = measurement[ialt, iact]['sun'][:]
            l2diag_gainCO2[ialt, iact, :] = l2product[ialt, iact]['gain_XCO2'][:]
            l2diag_gainH2O[ialt, iact, :] = l2product[ialt, iact]['gain_XH2O'][:]
            l2diag_gainCH4[ialt, iact, :] = l2product[ialt, iact]['gain_XCH4'][:]

    output_l2diag.close()

    return

def level1b_to_level2_processor(paths, global_config, local_config):

    # paths and local_config parameter

    l1b_path = paths.project + paths.data_interface + paths.interface_l1b
    afgl_path = paths.project + paths.data_afgl
    run_id = '_'+global_config['run_id']

    filename = local_config['filename']['level1b']+'_'+global_config['profile']+run_id

    # get the l1b files    
    l1b = get_l1b(l1b_path, filename)
    
    # get pixel mask
    if(local_config['pixel_mask']):
        print('take pixel mask')
        mask  = np.load(l1b_path+local_config['filename']['pixel_mask']+'.npy')
    else:       
        nact  = l1b['radiance'][0,:,0].size
        nwave = l1b['radiance'][0,0,:].size
        mask  = np.full((nact, nwave), True)
        
    # fig = plt.figure(figsize=(10, 8), dpi=100)
    # radiance = l1b['radiance'][0, :,:]
    # masked_data = np.ma.masked_where( np.invert(mask).data, radiance)
    # plt.imshow(masked_data, vmin = 1.E16,vmax = 4.E16, cmap = 'viridis', aspect = 4)

    # Internal lbl spectral grid
    wave_start = local_config['spec_settings']['wavestart']
    wave_end = local_config['spec_settings']['waveend']
    wave_extend = local_config['spec_settings']['wave_extend']
    dwave_lbl = local_config['spec_settings']['dwave']
    wave_lbl = np.arange(wave_start-wave_extend, wave_end+wave_extend, dwave_lbl)  # nm

    # Vertical layering
    nlay = local_config['atmosphere']['nlay']
    dzlay = local_config['atmosphere']['dzlay']
    psurf = local_config['atmosphere']['psurf']
    nlev = nlay + 1  # number of levels

    zlay = (np.arange(nlay-1, -1, -1)+0.5)*dzlay  # altitude of layer midpoint
    zlev = np.arange(nlev-1, -1, -1)*dzlay  # altitude of layer interfaces = levels

    # model atmosphere

    atm = libATM.atmosphere_data(zlay, zlev, psurf)
    atm.get_data_AFGL(paths.project + paths.data_afgl+local_config['std_atm'])

    # scale to some reference column mixing ratios
    xco2_ref = 405.  # ppm
    xco2 = np.sum(atm.CO2) / np.sum(atm.air) * 1.E6
    atm.CO2 = xco2_ref/xco2 * atm.CO2

    # Safe reference water, ch4 and co2 profiles
    ref_H2O = atm.H2O
    ref_CH4 = atm.CH4
    ref_CO2 = atm.CO2

    # XCO2 = np.sum(atm.CO2)/np.sum(atm.air)
    # XCH4 = np.sum(atm.CH4)/np.sum(atm.air)
    # XH2O = np.sum(atm.H2O)/np.sum(atm.air)
    # print(XCO2*1.E6, XCH4*1.E9, XH2O*1E6)

    sun = libRT.read_sun_spectrum_TSIS1HSRS(paths.project+paths.data_sol_spec+'/hybrid_reference_spectrum_c2021-03-04_with_unc.nc')
    sun_lbl = libRT.interpolate_sun(sun, wave_lbl)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(wave_lbl)
    molec.get_data_HITRAN('../data/hapi_data/', iso_ids)

    # Calculate optical properties
    xsec_file = paths.project + paths.data_tmp +  'optics_prop_l1bl2'+run_id+'.pkl'
    # If pickle file exists read from file
    if ((not os.path.exists(xsec_file)) or local_config['xsec_forced']):
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm)
        # Dump optics.prop dictionary into temporary pkl file
        pkl.dump(optics.prop, open(xsec_file, 'wb'))
    else:
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Read optics.prop dictionary from pickle file
        optics.prop = pkl.load(open(xsec_file, 'rb'))

    optics.set_opt_depth_species(atm, ['molec_01', 'molec_32', 'molec_07'])

    # Initialization of the least squares fit
    retrieval_init = {}
    retrieval_init['chi2 limit'] = local_config['retrieval_init']['chi2_lim']
    retrieval_init['maximum iteration'] = local_config['retrieval_init']['max_iter']
    retrieval_init['zlay'] = zlay
    retrieval_init['trace gases'] = {'CO2': {'init': 300,      'scaling': 1E-6, 'ref_profile': ref_CO2},
                                     'CH4': {'init': 1700,     'scaling': 1E-9, 'ref_profile': ref_CH4},
                                     'H2O': {'init': 1000,     'scaling': 1E-6, 'ref_profile': ref_H2O}}

#    retrieval_init['surface'] = {'alb0': 0.17}
    retrieval_init['wavelength lbl'] = wave_lbl
    retrieval_init['solar irradiance'] = sun_lbl

    nact = l1b['latitude'][0, :].size
    nalt = l1b['latitude'][:, 0].size

    l2product = {}
    xch4_model = 1.8E-6
    xco2_model = 405.E-6

    l2product = np.ndarray((nalt, nact), np.object_)
    measurement = np.ndarray((nalt, nact), np.object_)

    XCO2 = np.zeros(nact)

    for ialt in tqdm(range(nalt)):
        for iact in range(nact):
            wavelength = l1b['wavelength'][iact,mask[iact,:]].data
#            wavelength = l1b['wavelength'][iact, :]
            istart = np.argmin(np.abs(wavelength - wave_start))
            iend = np.argmin(np.abs(wavelength - wave_end))
            wave_meas = wavelength[istart:iend+1]  # nm
            
            nwave = wave_meas.size

            # define isrf object
            isrf = libNumTools.isrfct(wave_meas, wave_lbl)
            isrf.get_isrf(local_config['isrf_settings'])
           
            atm_ret = deepcopy(atm)  # to initialize each retrieval with the same atmosphere

            sun = isrf.isrf_convolution(sun_lbl)

            # Observation geometry
            measurement[ialt, iact] = {}
            measurement[ialt, iact]['wavelength'] = wave_meas[:]
            measurement[ialt, iact]['mu0'] = np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
            measurement[ialt, iact]['muv'] = np.cos(np.deg2rad(l1b['vza'][ialt, iact]))
            measurement[ialt, iact]['ymeas'] = l1b['radiance'][ialt, iact, mask[iact,:]][istart:iend+1].data
            measurement[ialt, iact]['Smeas'] = np.eye(nwave)*(l1b['noise'][ialt, iact, mask[iact,:]][istart:iend+1].data)**2
#            measurement[ialt, iact]['ymeas'] = l1b['radiance'][ialt, iact, istart:iend+1]
#            measurement[ialt, iact]['Smeas'] = np.eye(nwave)*(l1b['noise'][ialt, iact, istart:iend+1])**2
            measurement[ialt, iact]['sun'] = sun

            #derive first guess albedo from the maximum reflectance   
            ymeas_max=np.max(measurement[ialt, iact]['ymeas'])
            idx = np.where(measurement[ialt, iact]['ymeas'] == ymeas_max)[0][0]
            alb_first_guess = measurement[ialt, iact]['ymeas'][idx]/sun[idx]*np.pi/np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
            retrieval_init['surface'] = {'alb0': alb_first_guess, 'alb1': 0.0}

            # Non-scattering least squares fit
            l2product[ialt, iact] = libINV.Gauss_Newton_iteration(
                retrieval_init, atm_ret, optics, measurement[ialt,
                                                             iact], isrf, local_config['retrieval_init']['max_iter'],
                local_config['retrieval_init']['chi2_lim'])

            if(not l2product[ialt, iact]['convergence']):
                print('pixel did not converge (ialt,iact) = ', ialt,iact)
                
            # Define proxy product
            l2product[ialt, iact]['XCO2 proxy'] = l2product[ialt, iact]['XCO2']/l2product[ialt, iact]['XCH4']*xch4_model
            l2product[ialt, iact]['XCH4 proxy'] = l2product[ialt, iact]['XCH4']/l2product[ialt, iact]['XCO2']*xco2_model
            rel_error = np.sqrt((l2product[ialt, iact]['XCO2 precision']/l2product[ialt, iact]['XCO2'])**2 +
                                (l2product[ialt, iact]['XCH4 precision']/l2product[ialt, iact]['XCH4'])**2)
            l2product[ialt, iact]['XCO2 proxy precision'] = rel_error * l2product[ialt, iact]['XCO2']
            l2product[ialt, iact]['XCH4 proxy precision'] = rel_error * l2product[ialt, iact]['XCH4']

            XCO2[iact] = l2product[ialt, iact]['XCO2 proxy']
    # output to netcdf file

    print('=============================')
    print(np.mean(XCO2)*1.E6)
    print('=============================')
    
    l2_path = paths.project + paths.data_interface + paths.interface_l2 
    run_id = '_'+global_config['run_id']
    filename = local_config['filename']['level2'] + '_'+global_config['profile']+run_id
    level2_output(l2_path, filename, l2product, retrieval_init, l1b)
    # have to be fixed because of variable spectral size of measurement vector
#    filename = local_config['filename']['l2_output'] + '_'+global_config['profile']+'_diags_'+run_id
#    level2_diags_output(l2_path, filename, l2product, measurement)
    print('=> l1bl2 finished successfully for run_id ' + run_id + ' of profile '+ global_config['profile'] )

    return
