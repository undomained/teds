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
import multiprocessing
import datetime
import time
import h5py
import shutil
import logging
from itertools import repeat
import tqdm

from lib import libATM, libSGM, libRT_no2, libNumTools, constants
from lib.libWrite import writevariablefromname


class Dict2Class:
    """Convert a dictionaly to a class."""

    def __init__(self, arg_dict):
        self.__dict__.update(arg_dict)


def get_gm_data(filename):
    with nc.Dataset(filename, mode='r') as data:
        gm_data = {}
        gm_data['sza'] = deepcopy(data['sza'][:, :])
        gm_data['saa'] = deepcopy(data['saa'][:, :])
        gm_data['vza'] = deepcopy(data['vza'][:, :])
        gm_data['vaa'] = deepcopy(data['vaa'][:, :])
        gm_data['lat'] = deepcopy(data['lat'][:, :])
        gm_data['lon'] = deepcopy(data['lon'][:, :])
    return gm_data


def get_dem(fname, lat, lon, radius: int = 3):

    shp = lat.shape
    lats = lat.flatten()
    lons = lon.flatten()
    lons = np.where(lons<0., lons + 360., lons)

    with nc.Dataset(fname, 'r') as f:
        
        grp = f[f'DEM_RADIUS_{1000*radius:05d}']

        latitudes = grp['latdim'][:].data

        n_longitude = grp['n_longitude'][:].data
        n_longitude_sum = grp['n_longitude_sum'][:].data

        grid_idx = []
        sign = np.sign( latitudes[1] -  latitudes[0])
        lat_idx = np.round( np.interp(sign * lats, sign * latitudes, np.arange(latitudes.shape[0])) ).astype(int)
        for i in range(lats.shape[0]):
            longitude_idx = int(np.round(n_longitude[lat_idx[i]] * lons[i]/360.0))
            grid_idx.append( n_longitude_sum[lat_idx[i]] + longitude_idx )

        alt = grp['altitude'][:].data[grid_idx]

    alt[alt>1e19] = np.nan
    alt = np.reshape(alt, shp)

    return alt

def time_units_in_seconds(time_units):

    logger = logging.getLogger()


    tu = time_units.split(' ')[0]
    
    if tu == 'seconds':
        return 1.
    elif tu == 'minutes':
        return 60.
    elif tu == 'hours':
        return 3600.
    elif tu == 'days':
        return 24. * 3600.
    else:
        logger.error('Time unit not understood: ', tu)
        return(-1)

def get_cams_profiles( cfg, time, gm_data):

    logger = logging.getLogger()

    lats = gm_data['lat']
    lons = gm_data['lon']

    mmr_to_ppmv = {}
    mmr_to_ppmv['no2'] =  28.9644 / 46.0055 * 1e6
    mmr_to_ppmv['o3'] =  28.9644 / 47.9982 * 1e6

    cams = {}

    with nc.Dataset(cfg['atm']['cams']['path'], 'r') as f:
        time_units = f['time'].units
        epoch_date = np.array( (time_units.split(' ')[-2]).split('-')).astype(int)
        epoch_time = np.array( (time_units.split(' ')[-1]).split(':')).astype(int)
        epoch = datetime.datetime(epoch_date[0], epoch_date[1], epoch_date[2], epoch_time[0], epoch_time[1], epoch_time[2])

        d_time_hours = (time-epoch).total_seconds() / time_units_in_seconds(time_units) * np.ones(lats.shape)

        for key in f.variables.keys():
            cams[key] = f[key][:].data


    profile_data = {}
    profile_data['hyai'] = cams['hyai']
    profile_data['hybi'] = cams['hybi']
    profile_data['hyam'] = cams['hyam']
    profile_data['hybm'] = cams['hybm']

    # copy some parameters
   
    # if time only has one value, dont interpolate over time
    if len(cams['time']) == 1:
        time_slice = 0
        idx, w, shape = libNumTools.ndim_lin_interpol_get_indices_weights( [cams['lat'], cams['lon']], [lats, lons] )

    else:
        time_slice = slice(None)
        idx, w, shape = libNumTools.ndim_lin_interpol_get_indices_weights( [cams['time'], cams['lat'], cams['lon']], [d_time_hours, lats, lons] )

    # interpolate fields to time lat lon
    profile_data['psfc'] = libNumTools.ndim_lin_interpol_get_values( cams['SP'][time_slice,:,:], idx, w, shape) * 1e-2 # convert to hPa
    profile_data['psl'] = libNumTools.ndim_lin_interpol_get_values( cams['MSL'][time_slice,:,:], idx, w, shape) * 1e-2 # convert to hPa
    profile_data['zgeop'] = libNumTools.ndim_lin_interpol_get_values( cams['Z'][time_slice,:,:], idx, w, shape)  # hPa


    variables = cfg['atm']['gases'].copy()
    variables.append('t')

    for var in variables:
        if var == 'o3':
            varcams = 'go3'
        else:
            varcams = var

        if varcams not in cams:
            logger.error(f'Requested gas {var} not in CAMS file.')
            continue
        
        varlist = []
        for i in range( cams['t'].shape[1]):
            varlist.append( libNumTools.ndim_lin_interpol_get_values( cams[varcams][time_slice,i,:,:], idx, w, shape) )
        
        profile_data[var] = np.array( varlist )

        # convert mmr to ppmv
        if var in mmr_to_ppmv:
            profile_data[var] *= mmr_to_ppmv[var]

    return profile_data


def combine_mhh_cams( mhh_data, atm, species=['no2', 'co2', 'no']):

    # mhh in partial column, cams in ppmv
    # mhh at elevation grid, cams at hybrid sigma grid
    # mhh at higher resolution than cams
    # mhh in class, cams in dict
    # just add mhh regrid to cams --> assuming cams only contains background and microHH does not!!

    nalt, nact = atm['psfc'].shape

    # compute pressures for mhh elevations using 8km scale height and surface pressure
    mhh_plev = atm['psfc'][...,np.newaxis] * np.exp( -1 * mhh_data.zlev[np.newaxis,np.newaxis,...] / 8000.)
    mhh_play = atm['psfc'][...,np.newaxis] * np.exp( -1 * mhh_data.zlay[np.newaxis,np.newaxis,...] / 8000.)

    # convert mhh partial column [molec./m2] to ppmv

    # pressure drop per layer
    mhh_dplev = (mhh_plev[:,:,:-1] - mhh_plev[:,:,1:]) *1e2 # [Pa]

    for s in species:
        setattr( mhh_data, s + '_ppmv', getattr(mhh_data, s) * constants.g0 * constants.MDRYAIR * 1e6 / (constants.NA * mhh_dplev.data ) ) #[ppmv]

    # setattr( mhh_data, 'plev', mhh_plev )
    # setattr( mhh_data, 'play', mhh_play )

    # take mean over mhh layers inside cams layer

    nz = atm['hyai'].shape[0] # output number of layers

    for s in species:
        setattr( mhh_data, s + '_regridded', np.zeros((nz-1, nalt, nact)) )

    for ix in range(nalt):
        for iy in range(nact):
            preslev =  (1e-2 * atm['hyai'] + atm['psfc'][ix,iy] * atm['hybi'])

            for iz in range(nz-1):
                idx = np.where((mhh_play[ix,iy,:] >= preslev[iz]) & (mhh_play[ix,iy,:] < preslev[iz+1]))[0]

                if len(idx) < 1:
                    continue

                for s in species:
                    getattr(mhh_data, s+ '_regridded')[iz,ix,iy] = np.nanmean(getattr(mhh_data, s+'_ppmv')[ix,iy,idx])

    # combine mhh regrid and cams
     # just add mhh regrid to cams --> assuming cams only contains background and microHH does not!!
    for s in species:
        atm[s] += getattr(mhh_data, s+ '_regridded')

    return

def convert_atm_profiles(cfg, atm):
    # convert atmospheric profiles for usage in disamar

    profiles = {}

    if cfg['atm']['type'] == 'afgl':
         # afgl, input is in molec/m2, convert to ppmv
        nact, nalt=  atm.shape
        nlev =  atm[0, 0].plev.shape[0]

        variables = cfg['atm']['gases'].copy()
        variables.extend(['t','p'])
        for var in variables:
            profiles[var] = np.zeros((nact, nalt, nlev))

        for idx in range(nact):
            for idy in range(nalt):

                # psfc = atm[idx, idy].psurf / 100.0 # Pa --> hPa

                plev = atm[idx, idy].plev[::-1] / 100.0 # levels,  Pa --> hPa, surface --> TOA
                play = atm[idx, idy].play[::-1] / 100.0 # layers, Pa --> hPa, surface --> TOA

                tlev = atm[idx, idy].tlev[::-1]         # levels,  K --> hPa, surface --> TOA
                # tlay = atm[idx, idy].tlay[::-1]         # layers, K --> hPa, surface --> TOA

                air = atm[idx, idy].air[::-1]           # layers, molec/m2, surface --> TOA
                
                profiles['p'][idx, idy, :] = plev
                profiles['t'][idx, idy, :] = tlev

                for gas in cfg['atm']['gases']:
                    gas_prof = atm[idx, idy].NO2[::-1]           # layers, molec/m2, surface --> TOA

                    # convert molec/m2 to volume mixing ratio
                    gas_ppmv = gas_prof / air * 1e6 # molec/m2 --> ppmv

                    # interpolate layers to levels using log(p)

                    # no extrapolation, edge value is repeated
                    profiles[gas][idx, idy, :] = np.interp( np.log(plev[::-1]), np.log(play[::-1]), gas_ppmv[::-1] )[::-1] # strictly increasing

                    # extrapolation
                    # interp_gas = scipy.interpolate.RegularGridInterpolator(np.reshape(np.log(play.T),(1,len(play))), gas_ppmv, method='linear',bounds_error=False, fill_value=None)
                    # profiles[gas][idx, idy, :] = interp_gas(np.log(plev))

        # set last layer of profiles to TOA, otherwise disamar does not like it
        if (profiles['p'][:,:,-1] > 0.3 ).any():
            profiles['p'][:,:,-1] = 0.3
            profiles['t'][:,:,-1] = 250.0
            if 'no2' in cfg['atm']['gases']:
                profiles['no2'][:,:,-1] =6.2958549E-09
            if 'o3' in cfg['atm']['gases']:
                profiles['o3'][:,:,-1] = 7.4831730E-01


    elif cfg['atm']['type'] == 'cams':
        # cams, input is alread in molec/m2
        # create pressure profile and rearange arrays


        play = (1e-2 * atm['hyam'][np.newaxis,np.newaxis,::-1] + atm['psfc'][...,np.newaxis] * atm['hybm'][np.newaxis,np.newaxis,::-1]) # layers, hPa, surface --> TOA
        plev = (1e-2 * atm['hyai'][np.newaxis,np.newaxis,::-1] + atm['psfc'][...,np.newaxis] * atm['hybi'][np.newaxis,np.newaxis,::-1]) # levels, hPa, surface --> TOA

        # omit TOA layer CAMS, pressure is 0, disamar does not like it
        play = play[:,:,:-1]
        plev = plev[:,:,:-1]

        profiles['p'] = plev

        variables = cfg['atm']['gases'].copy()
        variables.append('t')

        nact, nalt, nlev =  plev.shape


        for var in variables:
            var_layer = np.transpose(atm[var][::-1,:,:], axes=[1,2,0]) # layers, ppmv, surface --> TOA

            # omit TOA layer CAMS
            var_layer = var_layer[:,:,:-1]

            profiles[var] = np.zeros((nact,nalt,nlev))
            for idx in range(nact):
                for idy in range(nalt):

                    # no extrapolation, edge value is repeated, get gas profile at levels
                    profiles[var][idx,idy,:] = np.interp( np.log(plev[idx,idy,::-1]), np.log(play[idx,idy,::-1]), var_layer[idx,idy,::-1] )[::-1] # strictly increasing

                    # extrapolation
                    # interp_var = scipy.interpolate.RegularGridInterpolator(np.reshape(np.log(play[idx,idy,:].T),(1,len(play[idx,idy,:]))), var_layer[idx,idy,:], method='linear',bounds_error=False, fill_value=None)
                    # profiles[var][idx,idy,:] = interp_gas(np.log(plev[idx,idy,:]))

    return profiles



def set_disamar_cfg_sim(cfg, dis_cfg, ground_points, profiles, albedo, i_t, i_x):

    # adapted from Pepijn's E2E
    # Modify the disamar input file for simulation of spectra
    # all instruments modification off, no slit, high resolution

    logger = logging.getLogger()

    # GENERAL

    if cfg['rtm']['dismas_sim']:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(1)
        dis_cfg['GENERAL','method', 'ignoreSlitSim'].setvalue(0)
        dis_cfg['INSTRUMENT','slit_parameters', 'FWHM_irradiance_sim'].setvalue(cfg['rtm']['dismas_fwhm'])
        dis_cfg['INSTRUMENT','slit_parameters', 'FWHM_radiance_sim'].setvalue(cfg['rtm']['dismas_fwhm'])

    else:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(0)
        dis_cfg['GENERAL','method', 'ignoreSlitSim'].setvalue(1)


    # get datetime
    # dt = ground_points['epoch'] + datetime.timedelta( seconds=ground_points['seconds_from_epoch'][i_t] )

    # add the lat-lon information
    dis_cfg['GENERAL','external_data', 'latitude'].setvalue([ground_points['lat'][i_t,i_x]])
    dis_cfg['GENERAL','external_data', 'longitude'].setvalue([ground_points['lon'][i_t,i_x]])
    dis_cfg['GENERAL','external_data', 'AtrackNumber'].setvalue([i_t])
    dis_cfg['GENERAL','external_data', 'XtrackNumber'].setvalue([i_x])
    # dis_cfg['GENERAL','external_data', 'year'].setvalue([dt.year])
    # dis_cfg['GENERAL','external_data', 'month'].setvalue([dt.month])
    # dis_cfg['GENERAL','external_data', 'day'].setvalue([dt.day])
    # dis_cfg['GENERAL','external_data', 'hour'].setvalue([dt.hour])
    # dis_cfg['GENERAL','external_data', 'second'].setvalue([dt.second])

    # INSTRUMENT
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_start'].setvalue(cfg['rtm']['wave_start'])
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_end'].setvalue(cfg['rtm']['wave_end'])
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_step'].setvalue(cfg['rtm']['dwave'])

    # GEOMETRY   
    dis_cfg['GEOMETRY','geometry', 'solar_zenith_angle_sim'].setvalue( ground_points['sza'][i_t,i_x] )
    dis_cfg['GEOMETRY','geometry', 'instrument_nadir_angle_sim'].setvalue( ground_points['vza'][i_t,i_x] )
    dis_cfg['GEOMETRY','geometry', 'solar_azimuth_angle_sim'].setvalue( ground_points['saa'][i_t,i_x] )
    dis_cfg['GEOMETRY','geometry', 'instrument_azimuth_angle_sim'].setvalue( ground_points['vaa'][i_t,i_x] )
    

    # Profiles
    # PT profile
    pt_sim = dis_cfg['PRESSURE_TEMPERATURE', 'PT_sim', 'PT']
    pt_sim.set_rawvalue(np.asarray([profiles['p'][i_t,i_x,:], profiles['t'][i_t,i_x,:]]).T)
    pt_retr = dis_cfg['PRESSURE_TEMPERATURE', 'PT_retr', 'PT']
    pt_retr.set_rawvalue(np.asarray([profiles['p'][i_t,i_x,:], profiles['t'][i_t,i_x,:], np.ones(profiles['p'][i_t,i_x,:].shape, dtype=float)]).T)

    # gas profiles
    for gas in cfg['atm']['gases']:
        gas_vmr = np.asarray([profiles['p'][i_t,i_x,:], profiles[gas][i_t,i_x,:]]).T
        gas_vmr_error = np.asarray([profiles['p'][i_t,i_x,:], profiles[gas][i_t,i_x,:], np.ones(profiles['p'][i_t,i_x,:].shape, dtype=float) * 20.]).T
        dis_cfg[gas.upper(), 'profile', 'P_vmr_ppmv_sim'].set_rawvalue(gas_vmr)
        dis_cfg[gas.upper(), 'profile', 'P_vmr_ppmv_error_percent_retr'].set_rawvalue(gas_vmr_error)
        
    # SURFACE

    # find albedo wvl
    match cfg['S2_albedo']['band']:
        case 'B01':
            albedo_wvl = 442.5 # [nm]
        case 'B02':
            albedo_wvl = 492.3
        case _:
            logger.error('unknown S2 albedo band')
            albedo_wvl = 440.0

    dis_cfg['SURFACE', 'wavelDependentSim', 'wavelSurfAlbedo'].setvalue(albedo_wvl)
    dis_cfg['SURFACE', 'wavelDependentSim', 'surfAlbedo'].setvalue(albedo[i_t,i_x])

    dis_cfg['SURFACE','pressure', 'surfPressureSim'].setvalue(profiles['p'][i_t,i_x, 0])
    dis_cfg['SURFACE','pressure', 'surfPressureRetr'].setvalue(profiles['p'][i_t,i_x, 0])


    return dis_cfg



def run_disamar(filename,disamar_exe):


    dis_cfg = libRT_no2.RT_configuration(filename=filename)
    output_filename = filename.replace('.in', '.h5')

    RT = libRT_no2.rt_run(cfg=dis_cfg, disamar=disamar_exe, output=output_filename, quiet=True, debug=False)

    try:
        starttime = time.time()
        RT()
        # logging.info(f'finished: {filename} in {np.round(time.time()-starttime,1)} s')
        return 0
    except:
        # logging.error(f'failed: {filename}')
        return -1


def read_disamar_output(gm_data,tmp_dir):

    logger = logging.getLogger()

    nalt, nact = gm_data['lat'].shape

    dis_output = {}

    for iact in range(nact):
        for ialt in range(nalt):
            file = f'{tmp_dir}/alt{ialt:04d}_act{iact:03d}.h5'
            
            if os.path.isfile(file) is False:
                logger.error(f'Disamar output file {file} not found')
                continue

            with h5py.File(file, 'r') as f:

                if 'wavelength_lbl' not in dis_output:
                    dis_output['wavelength_lbl'] = f['specifications_sim/wavelength_radiance_band_1'][:]            # nm
                    dis_output['solar_irradiance'] = f['radiance_and_irradiance/solar_irradiance_band_1'][:] *1.e4  # ph/s/nm/cm2 --> ph/s/nm/m2/sr
                    
                    nwvl = len(dis_output['wavelength_lbl'])
                    dis_output['radiance'] = np.full((nalt,nact,nwvl), np.nan)

                dis_output['radiance'][ialt,iact,:] = f['radiance_and_irradiance/earth_radiance_band_1'][:] *1.e4   # ph/s/nm/cm2/sr --> ph/s/nm/m2/sr


    if 'wavelength_lbl' not in dis_output:
        logger.error('No DISAMAR output found. Exiting')
        sys.exit()

    return dis_output


def sgm_output_radio(config, rad_output):
    # write radiances
    nalt, nact, nlbl = rad_output['radiance'].shape
    # open file
    output_rad = nc.Dataset(config['sgm_rad_file'], mode='w')
    output_rad.title = 'Tango E2ES SGM radiometric scene'
    # output_rad.config = str(config)
    output_rad.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')

    output_rad.createDimension('bins_spectral', nlbl)     # spectral axis
    output_rad.createDimension('bins_across_track', nact)     # across track axis
    output_rad.createDimension('bins_along_track', nalt)     # along track axis
    # wavelength
    _ = writevariablefromname(output_rad, 'wavelength', ('bins_spectral',), rad_output['wavelength_lbl'])
    # solar irradiance
    _ = writevariablefromname(output_rad, 'solarirradiance', ('bins_spectral',), rad_output['solar_irradiance'])
    # radiance
    dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    _ = writevariablefromname(output_rad, 'radiance_sgm', dims, rad_output['radiance'])
    output_rad.close()


def sgm_output_atm_afgl(filename_atm, atm, albedo, gm_data, microhh_data):
    # write afgl atm to netcdf
    
    gases = config['atm']['gases']

    nalt, nact = gm_data["lat"].shape

    nlay, nlev = atm[0,0].zlay.size, atm[0,0].zlev.size

    # write 3d data
    # write the variables to an array
    zlev = np.zeros((nalt, nact, nlev))
    zlay = np.zeros((nalt, nact, nlay))
    dcolco2 = np.zeros((nalt, nact, nlay))
    dcolch4 = np.zeros((nalt, nact, nlay))
    dcolh2o = np.zeros((nalt, nact, nlay))
    dcolno2 = np.zeros((nalt, nact, nlay))
    dcolo3 = np.zeros((nalt, nact, nlay))
    ppmv_no2 = np.zeros((nalt, nact, nlay))
    ppmv_o3 = np.zeros((nalt, nact, nlay))

    # Total column integrated values of CO2, CH4, H2O, NO2 and air
    xco2 = np.zeros((nalt, nact))
    xch4 = np.zeros((nalt, nact))
    xh2o = np.zeros((nalt, nact))
    col_air = np.zeros((nalt, nact))
    col_no2 = np.zeros((nalt, nact))
    col_o3 = np.zeros((nalt, nact))

    for ialt in range(nalt):
        for iact in range(nact):
            zlay[ialt, iact, :] = atm[ialt][iact].zlay
            zlev[ialt, iact, :] = atm[ialt][iact].zlev
            dcolco2[ialt, iact, :] = atm[ialt][iact].CO2[:] # [molec/m2]
            dcolch4[ialt, iact, :] = atm[ialt][iact].CH4[:] # [molec/m2]
            dcolh2o[ialt, iact, :] = atm[ialt][iact].H2O[:] # [molec/m2]
            dcolno2[ialt, iact, :] = atm[ialt][iact].NO2[:]*1e-4 # [molec/cm2]
            dcolo3[ialt, iact, :] = atm[ialt][iact].O3[:]*1e-4 # [molec/cm2]

            ppmv_no2[ialt,iact,:] = atm[ialt][iact].NO2[:]/atm[ialt, iact].air[:]*1.e6 # [ppmv]
            ppmv_o3[ialt,iact,:] = atm[ialt][iact].O3[:]/atm[ialt, iact].air[:]*1.e6 # [ppmv]

            XAIR = np.sum(atm[ialt, iact].air[:])
            xco2[ialt, iact] = np.sum(atm[ialt, iact].CO2[:])/XAIR*1.e6  # [ppmv]
            xch4[ialt, iact] = np.sum(atm[ialt, iact].CH4[:])/XAIR*1.e9  # [ppbv]
            xh2o[ialt, iact] = np.sum(atm[ialt, iact].H2O[:])/XAIR*1.e6  # [ppmv]
            col_no2[ialt, iact] = np.sum(atm[ialt, iact].NO2[:])*1e-4    # [molec/cm2]
            col_o3[ialt, iact] = np.sum(atm[ialt, iact].O3[:])*1e-4    # [molec/cm2]
            col_air[ialt, iact] = XAIR

    # init file
    output_atm = nc.Dataset(filename_atm, mode='w')
    output_atm.title = 'Tango E2ES SGM atmospheric scene'
    # output_atm.config = str(config)
    output_atm.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')

    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis


    dims_lay = ('bins_along_track', 'bins_across_track', 'number_layers')
    dims_lev = ('bins_along_track', 'bins_across_track', 'number_levels')
    dims_2d = ('bins_along_track', 'bins_across_track')


    if ("co2" in gases):
        # columndensity_co2
        _ = writevariablefromname(output_atm, 'subcol_density_co2', dims_lay, dcolco2)
        # column_co2
        var_co2 = writevariablefromname(output_atm, 'XCO2', dims_2d, xco2)
    if ("ch4" in gases):
        # columndensity_ch4
        _ = writevariablefromname(output_atm, 'subcol_density_ch4', dims_lay, dcolch4)
        # column_ch4
        var_ch4 = writevariablefromname(output_atm, 'XCH4', dims_2d, xch4)
    if ("ch4" in gases) or ("co2" in gases) :
        # columndensity_h2o
        _ = writevariablefromname(output_atm, 'subcol_density_h2o', dims_lay, dcolh2o)
    if ("no2" in gases):
        # columndensity_no2
        _ = writevariablefromname(output_atm, 'subcol_density_no2', dims_lay, dcolno2)
        # concentration_no2
        _ = writevariablefromname(output_atm, 'conc_no2', dims_lay, ppmv_no2)
        # total column no2
        var_no2 = writevariablefromname(output_atm, 'column_no2', dims_2d, col_no2)
    if ("o3" in gases):
        # columndensity_o3
        _ = writevariablefromname(output_atm, 'subcol_density_o3', dims_lay, dcolo3)
        _ = writevariablefromname(output_atm, 'conc_o3', dims_lay, ppmv_o3)
        _ = writevariablefromname(output_atm, 'column_o3', dims_2d, col_o3)
    if ("ch4" in gases) or ("co2" in gases) :
        # column_h2o
        _ = writevariablefromname(output_atm, 'XH2O', dims_2d, xh2o)
        # column_air
        _ = writevariablefromname(output_atm, 'column_air', dims_2d, col_air)

    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', dims_lay, zlay)
    # level height
    _ = writevariablefromname(output_atm, 'levelheight', dims_lev, zlev)
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', dims_2d, albedo)
    # lat lon
    _ = writevariablefromname(output_atm, 'latitude', dims_2d, gm_data["lat"])
    _ = writevariablefromname(output_atm, 'longitude', dims_2d, gm_data["lon"])


    # write microhh attributes
    if microhh_data is not None:
        if ("co2" in gases):
            var_co2.setncattr("source", microhh_data.__getattribute__("co2_source"))
            var_co2.setncattr("emission_kgps", microhh_data.__getattribute__("co2_emission_in_kgps"))

        if "ch4" in gases:
            var_ch4.setncattr("source", microhh_data.__getattribute__("ch4_source"))
            var_ch4.setncattr("emission_kgps", microhh_data.__getattribute__("ch4_emission_in_kgps"))

        # column_no2
        if "no2" in gases:
            var_no2.setncattr("source", microhh_data.__getattribute__("no2_source"))
            var_no2.setncattr("emission_kgps", microhh_data.__getattribute__("no2_emission_in_kgps"))

    output_atm.close()

    return


def sgm_output_atm_cams(config, atm, albedo, gm_data, microhh_data):
    # write cams atm to netcdf

    gases = config['atm']['gases']

    nalt, nact = gm_data["lat"].shape
    nlay, nlev = len(atm['hyam']), len(atm['hyai'])

    # calculate pressure profiles

    play = (1e-2 * atm['hyam'][np.newaxis,np.newaxis,...] + atm['psfc'][...,np.newaxis] * atm['hybm'][np.newaxis,np.newaxis,...]) # layers, hPa, TOA --> surface
    plev = (1e-2 * atm['hyai'][np.newaxis,np.newaxis,...] + atm['psfc'][...,np.newaxis] * atm['hybi'][np.newaxis,np.newaxis,...]) # levels, hPa, TOA --> surface

    # create file
    output_atm = nc.Dataset(config['sgm_atm_file'], mode='w')
    output_atm.title = 'Tango E2ES SGM atmospheric scene'
    # output_atm.config = str(config)
    output_atm.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')


    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    # add variables
    dims_layer = ('bins_along_track', 'bins_across_track', 'number_layers')
    dims_level = ('bins_along_track', 'bins_across_track', 'number_levels')
    dims_2d = ('bins_along_track', 'bins_across_track')

    _ = writevariablefromname(output_atm, 'latitude', dims_2d, gm_data["lat"])
    _ = writevariablefromname(output_atm, 'longitude', dims_2d, gm_data["lon"])
    _ = writevariablefromname(output_atm, 'pressure_levels', dims_level, plev)
    _ = writevariablefromname(output_atm, 'pressure_layers', dims_layer, play)
    _ = writevariablefromname(output_atm, 'temperature', dims_layer, np.transpose(atm['t'], axes=[1,2,0]))
    _ = writevariablefromname(output_atm, 'albedo', dims_2d, albedo)

    # pressure drop per layer [Pa]
    pdlev = (plev[:,:,1:] - plev[:,:,:-1])*1e2

    for gas in gases:

        # ppmv profile
        gas_ppmv = np.transpose(atm[gas], axes=[1,2,0]) # [ppmv]
        _ = writevariablefromname(output_atm, 'conc_'+ gas, dims_layer, gas_ppmv) # [ppmv]

        # Convert volume mixing ratio to NO2 partial column

        # partial column in [mol/m2] at layers
        gas_partialcolumn = pdlev * gas_ppmv*1e-6 / ( constants.g0 * constants.MDRYAIR ) 
        # convert to [molec/cm^2]
        gas_partialcolumn *=  constants.NA * 1e-4

        # calculate total column
        gas_totalcolumn = gas_partialcolumn.sum(axis=-1)

        # write partial and total column
        var_subcol = writevariablefromname(output_atm, 'subcol_density_'+gas, dims_layer, gas_partialcolumn) # [molec/cm2]
        var_totalcolumn = writevariablefromname(output_atm, 'column_'+gas, dims_2d, gas_totalcolumn) # [molec/cm2]

        # write microHH attributes
        if config['atm']['microHH']['use']:
            if gas not in config['atm']['microHH']['gases']:
                continue
            var_totalcolumn.setncattr("source", microhh_data.__getattribute__(f"{gas}_source"))
            var_totalcolumn.setncattr("emission_kgps", microhh_data.__getattribute__(f"{gas}_emission_in_kgps"))


    output_atm.close()

    return




def scene_generation_module_nitro(logging, config):
    """
    Scene generation algorithm for NO2

    Note: for now only works for profile: orbit
    """
    logging.info('Starting SGM calculation')

    start_time = time.time()

    # =============================================================================================
    # 1) get the geometry data
    # =============================================================================================

    gm_data = get_gm_data(config['gm_file'])
    nalt, nact = gm_data['sza'].shape
    logging.info(f"Pixels along track: {nalt} ; Pixels across track: {nact}")

    # =============================================================================================
    # 2) get collocated S2 albedo data
    # =============================================================================================

    albedo = np.zeros([nalt, nact])

    file_exists = os.path.isfile(config['S2_albedo']['dump'])
    if (file_exists and (not config['S2_albedo']['forced'])):
        logging.info(f"Loading S2 data from dump file: {config['S2_albedo']['dump']}")
        albedo = np.load(config['S2_albedo']['dump'])
    else:
        logging.info(f"Downloading S2 data for band: {config['S2_albedo']['band']}")
        albedo = libSGM.get_sentinel2_albedo(gm_data['lat'], gm_data['lon'], config,band=config['S2_albedo']['band'])
        np.save(config['S2_albedo']['dump'], albedo)    # .npy extension is added if not given
    
    # clip albedo
    albedo[albedo>1.0] = 1.0
    albedo[albedo<0.0] = 0.0

    # =============================================================================================
    # 3) get collocated microHH data and apply convolution
    # =============================================================================================
    if config['atm']['microHH']['use']:
        if ((not os.path.exists(config['atm']['microHH']['dump'])) or config['atm']['microHH']['forced']):
            logging.info(f"Loading microHH data: {config['atm']['microHH']['path_data']}")

            microhh_data = libATM.get_atmosphericdata(gm_data['lat'], gm_data['lon'], config['atm']['microHH'],
                                                   config['kernel_parameter'])
            # Dump microHH dictionary into temporary pkl file
            pickle.dump(microhh_data.__dict__, open(config['atm']['microHH']['dump'], 'wb'))
        else:
            logging.info(f"Loading microHH data from dump file: {config['atm']['microHH']['dump']}")
            # Read microHH from pickle file
            microhh_data = Dict2Class(pickle.load(open(config['atm']['microHH']['dump'], 'rb')))
    else:
        microhh_data = None

    # =============================================================================================
    # 4) build model atmosphere and write output
    # =============================================================================================
    match config['atm']['type']:

        case 'afgl':

            logging.info(f"Loading afgl atmosphere: {config['atm']['afgl']['path']}")

            # 4A) get a model atm from AFGL files

            nlay = config['atm']['afgl']['nlay']  # number of layers
            dzlay = config['atm']['afgl']['dzlay']
            # we assume the same standard atm for all pixels of the granule
            atm = libATM.get_AFGL_atm_homogenous_distribution(config['atm']['afgl']['path'], nlay, dzlay)
            
            if config['atm']['microHH']['use']:
                #combine the microHH data with standard atm afgl, regridding microHH to std atm layers
                atm = libATM.combine_meteo_standard_atm(microhh_data, atm, config['atm']['microHH']['gases'])
            else:
                atm = np.full((nalt,nact), atm, dtype=object)

            sgm_output_atm_afgl(config['sgm_atm_file'], atm, albedo, gm_data, microhh_data)
            
        case 'cams':

            logging.info(f"Loading CAMS atmosphere: {config['atm']['cams']['path']}")


            # 4B) get atm from CAMS

            # interpolate CAMS field
            atm = get_cams_profiles(config, config['atm']['cams']['start'], gm_data)

            if config['atm']['dem']['use']:

                logging.info(f"Loading DEM: {config['atm']['dem']['path']}")


                # correct surface pressure with DEM
                atm['zsfc'] = get_dem(config['atm']['dem']['path'], gm_data['lat'], gm_data['lon'])
                atm['psfc'] = atm['psl'] * np.exp( -1 * atm['zsfc'] / 8000.) # use 8 km scale height

            if config['atm']['microHH']['use']:
                # combine the microHH data with CAMS atm, regridding microHH to CAMS layers
                combine_mhh_cams( microhh_data, atm, species=['no2'])


            logging.info(f"Writing scene atmosphere: {config['sgm_atm_file']}")


            sgm_output_atm_cams(config, atm, albedo, gm_data, microhh_data)


    # =============================================================================================
    # 5) radiative transfer simulations with DISAMAR
    # =============================================================================================

    logging.info('Radiative transfer simulation')

    # 5A) generate the disamar config files and tmp dirs
    logging.info('Creating config files DISAMAR')

    # convert atm profiles to disamar format
    dis_profiles = convert_atm_profiles(config, atm)

    if os.path.isfile(config['rtm']['disamar_cfg_template']):
        dis_cfg = libRT_no2.RT_configuration(filename=config['rtm']['disamar_cfg_template'])
    else:
        logging.error(f"File {config['rtm']['disamar_cfg_template']} not found")

    dis_cfg_filenames=[]

    timestamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    tmp_dir = '{}/{}'.format( config['rtm']['tmp_dir'], timestamp)
    if not os.path.isdir(tmp_dir):
        logging.info(f'Creating tmp directory DISAMAR: {tmp_dir}')

        os.makedirs(tmp_dir, exist_ok=True)

    for ialt in range(nalt):
        for iact in range(nact):

            dis_cfg = set_disamar_cfg_sim(config, dis_cfg, gm_data, dis_profiles, albedo, ialt, iact)

            filename = f'{tmp_dir}/alt{ialt:04d}_act{iact:03d}.in'
            dis_cfg.write(filename=filename)
            dis_cfg_filenames.append(filename)
    

    # 5B) run disamar in parallel
    logging.info('Running DISAMAR')

    with multiprocessing.Pool(config['rtm']['n_threads']) as pool:
        # stat = pool.starmap(run_disamar, zip(dis_cfg_filenames, repeat(config['rtm']['disamar_exe'])),chunksize=1)
        stat = pool.starmap(run_disamar, tqdm.tqdm(zip(dis_cfg_filenames, repeat(config['rtm']['disamar_exe'])),total=len(dis_cfg_filenames)),chunksize=1)

    if sum(stat) > 0:
        logging.error('Error in at least one RT calculation run with DISAMAR')

    # 5C) read disamar output
    logging.info('Reading DISAMAR output')

    dis_output = read_disamar_output(gm_data, tmp_dir)

    # cleanup
    if config['rtm']['cleanup']:
        logging.info('Cleaning up tmp directory DISAMAR')
        shutil.rmtree(tmp_dir)

    # =============================================================================================
    # 6) sgm output to radiometric file
    # =============================================================================================
    logging.info(f"Writing scene radiance to: {config['sgm_rad_file']}")

    sgm_output_radio(config, dis_output)

    logging.info(f'SGM calculation finished in {np.round(time.time()-start_time,1)} s')
    return


if __name__ == '__main__':
    
    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )


    # call with:
    # python sgm_no2.py sgm.yaml

    # or with logging to file:
    # python sgm_no2.py sgm.yaml sgm.log


    # setup the logging to screen and to file
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

    if len(sys.argv) > 2:
        fh = logging.FileHandler(sys.argv[2], mode='w')
        fh.setLevel(logging.ERROR)
        fh.setFormatter(formatter)

        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter) 

        logging.basicConfig(level=logging.INFO, handlers = [ch,fh])

        logging.info(f'Logging to file: {sys.argv[2]}')

    else:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter) 
        logging.basicConfig(level=logging.INFO,handlers = [ch])

    logging.info(f'Reading config file: {sys.argv[1]}')
    config = yaml.safe_load(open(sys.argv[1]))
    scene_generation_module_nitro(logging,config)

