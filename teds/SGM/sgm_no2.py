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
from tqdm import tqdm
from lib.libWrite import writevariablefromname
import scipy.interpolate
import multiprocessing
import datetime

from lib import libATM, libSGM, libRT_no2, libNumTools


class Dict2Class:
    """Convert a dictionary to a class."""

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


def sgm_output_atm_afgl(filename_atm, atm, albedo, gm_data, meteodata=None, gases=None):
    # write afgl atm to netcdf


    nalt, nact = gm_data["lat"].shape
    # write atmosphere
    nlay, nlev = atm[0][0].zlay.size, atm[0][0].zlev.size
    # file
    output_atm = nc.Dataset(filename_atm, mode='w')
    output_atm.title = 'Tango E2ES SGM atmospheric scene'
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
    dcolno2 = np.zeros((nalt, nact, nlay))
    ppmv_no2 = np.zeros((nalt, nact, nlay))

    for ialt in range(nalt):
        for iact in range(nact):
            zlay[ialt, iact, :] = atm[ialt][iact].zlay
            zlev[ialt, iact, :] = atm[ialt][iact].zlev
            dcolco2[ialt, iact, :] = atm[ialt][iact].CO2[:]
            dcolch4[ialt, iact, :] = atm[ialt][iact].CH4[:]
            dcolh2o[ialt, iact, :] = atm[ialt][iact].H2O[:]
            dcolno2[ialt, iact, :] = atm[ialt][iact].NO2[:]

            ppmv_no2[ialt,iact,:] = atm[ialt][iact].NO2[:]/atm[ialt, iact].air[:]*1.e6


    _dims = ('bins_along_track', 'bins_across_track', 'number_layers')
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims, zlay)
    if ("co2" in gases):
        # columndensity_co2
        _ = writevariablefromname(output_atm, 'subcol_density_co2', _dims, dcolco2)
    if ("ch4" in gases):
        # columndensity_ch4
        _ = writevariablefromname(output_atm, 'subcol_density_ch4', _dims, dcolch4)
    if ("ch4" in gases) or ("co2" in gases) :
        # columndensity_h2o
        _ = writevariablefromname(output_atm, 'subcol_density_h2o', _dims, dcolh2o)
    if ("no2" in gases):
        # columndensity_no2
        _ = writevariablefromname(output_atm, 'subcol_density_no2', _dims, dcolno2*1e-4) # [molec/cm2]
        _ = writevariablefromname(output_atm, 'no2_profile', _dims, ppmv_no2) # [ppmv]

    # level height
    _dims = ('bins_along_track', 'bins_across_track', 'number_levels')
    _ = writevariablefromname(output_atm, 'levelheight', _dims, zlev)

    # Total column integrated values of CO2, CH4, H2O, NO2 and air
    xco2 = np.zeros((nalt, nact))
    xch4 = np.zeros((nalt, nact))
    xh2o = np.zeros((nalt, nact))
    col_air = np.zeros((nalt, nact))
    col_no2 = np.zeros((nalt, nact))
    for ialt in range(nalt):
        for iact in range(nact):
            XAIR = np.sum(atm[ialt, iact].air[:])
            xco2[ialt, iact] = np.sum(atm[ialt, iact].CO2[:])/XAIR*1.e6  # [ppmv]
            xch4[ialt, iact] = np.sum(atm[ialt, iact].CH4[:])/XAIR*1.e9  # [ppbv]
            xh2o[ialt, iact] = np.sum(atm[ialt, iact].H2O[:])/XAIR*1.e6  # [ppmv]
            col_no2[ialt, iact] = np.sum(atm[ialt, iact].NO2[:])*1e-4    # [molec/cm2]
            col_air[ialt, iact] = XAIR
    _dims = ('bins_along_track', 'bins_across_track')
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', _dims, albedo)

    # write new attributes
    if gases is not None:
        if ("co2" in gases):
            # column_co2
            var_co2 = writevariablefromname(output_atm, 'XCO2', _dims, xco2)
            var_co2.setncattr("source", meteodata.__getattribute__("co2_source"))
            var_co2.setncattr("emission_kgps", meteodata.__getattribute__("co2_emission_in_kgps"))
    # column_ch4
    if (gases is not None):
        if "ch4" in gases:
            var_ch4 = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
            var_ch4.setncattr("source", meteodata.__getattribute__("ch4_source"))
            var_ch4.setncattr("emission_kgps", meteodata.__getattribute__("ch4_emission_in_kgps"))

    # column_no2
    if (gases is not None):
        if "no2" in gases:
            var_no2 = writevariablefromname(output_atm, 'column_no2', _dims, col_no2)
            var_no2.setncattr("source", meteodata.__getattribute__("no2_source"))
            var_no2.setncattr("emission_kgps", meteodata.__getattribute__("no2_emission_in_kgps"))

    if ("ch4" in gases) or ("co2" in gases) :
        # column_h2o
        _ = writevariablefromname(output_atm, 'XH2O', _dims, xh2o)
        # column_air
        _ = writevariablefromname(output_atm, 'column_air', _dims, col_air)

    # add coordinates to SGM atmosphere
    _ = writevariablefromname(output_atm, 'latitude', _dims, gm_data["lat"])
    _ = writevariablefromname(output_atm, 'longitude', _dims, gm_data["lon"])
    output_atm.close()

    return


def sgm_output_atm_cams(config, atm, albedo, gm_data):
    # write cams atm to netcdf

    nalt, nact = gm_data["lat"].shape
    nlay, nlev = len(atm['hyam']), len(atm['hyai'])

    # calculate pressure profiles

    play = (1e-2 * atm['hyam'][np.newaxis,np.newaxis,::-1] + atm['psfc'][...,np.newaxis] * atm['hybm'][np.newaxis,np.newaxis,::-1]) # layers, hPa, surface --> TOA
    plev = (1e-2 * atm['hyai'][np.newaxis,np.newaxis,::-1] + atm['psfc'][...,np.newaxis] * atm['hybi'][np.newaxis,np.newaxis,::-1]) # levels, hPa, surface --> TOA

    # create file
    output_atm = nc.Dataset(config['output']['atm'], mode='w')
    output_atm.title = 'Tango E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    # add variables
    _dims_layer = ('bins_along_track', 'bins_across_track', 'number_layers')
    _dims_level = ('bins_along_track', 'bins_across_track', 'number_levels')
    _dims_surface = ('bins_along_track', 'bins_across_track')

    _ = writevariablefromname(output_atm, 'latitude', _dims_surface, gm_data["lat"])
    _ = writevariablefromname(output_atm, 'longitude', _dims_surface, gm_data["lon"])
    _ = writevariablefromname(output_atm, 'pressure_levels', _dims_level, plev)
    _ = writevariablefromname(output_atm, 'pressure_layers', _dims_layer, play)
    _ = writevariablefromname(output_atm, 'temperature', _dims_layer, np.transpose(atm['t'][::-1], axes=[1,2,0]))
    _ = writevariablefromname(output_atm, 'albedo', _dims_surface, albedo)

    for gas in config['atm']['gases']:
        _ = writevariablefromname(output_atm, gas+'_profile', _dims_layer, np.transpose(atm[gas][::-1], axes=[1,2,0])) # [ppmv]

# >>>>>>>>>>> to do here:  calculate total column




    output_atm.close()

    return


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
        print('Time unit not understood: ', tu)
        return(-1)

def get_cams_profiles( cfg, time, gm_data):

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

    # interpolate fields to time lat lon
    idx, w, shape = libNumTools.ndim_lin_interpol_get_indices_weights( [cams['time'], cams['lat'], cams['lon']], [d_time_hours, lats, lons] )
    profile_data['psfc'] = libNumTools.ndim_lin_interpol_get_values( cams['SP'][:,:,:], idx, w, shape) * 1e-2 # convert to hPa
    profile_data['psl'] = libNumTools.ndim_lin_interpol_get_values( cams['MSL'][:,:,:], idx, w, shape) * 1e-2 # convert to hPa
    profile_data['zgeop'] = libNumTools.ndim_lin_interpol_get_values( cams['Z'][:,:,:], idx, w, shape)  # convert to hPa

    # interpolate fields to lat lon
    idx, w, shape = libNumTools.ndim_lin_interpol_get_indices_weights( [cams['lat'], cams['lon']], [lats, lons] )

    variables = cfg['atm']['gases'].copy()
    variables.append('t')

    for var in variables:
        if var not in cams:
            print(f'Error, requested gas {var} not in CAMS file.')
            continue
        
        varlist = []
        for i in range( cams['t'].shape[0]):
            varlist.append( libNumTools.ndim_lin_interpol_get_values( cams[var][i,:,:], idx, w, shape) )
        
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

    # constants
    M_air   = 28.9644e-3                    # molar mass air [kg/mol]
    g       = 9.80665                       # gravitational constant Earth [kg/s^2]
    N_A     = 6.02214076e23                 # Avogadros number [mol-1]

    # pressure drop per layer
    mhh_dplev = (mhh_plev[:,:,:-1] - mhh_plev[:,:,1:]) *1e2 # [Pa]

    for s in species:
        setattr( mhh_data, s + '_ppmv', getattr(mhh_data, s) * g * M_air * 1e6 / (N_A * mhh_dplev.data ) ) #[ppmv]

    # setattr( mhh_data, 'plev', mhh_plev )
    # setattr( mhh_data, 'play', mhh_play )

    # take mean over mhh layers inside cams layer

    nz = atm['hyai'].shape[0] # output number of layers

    d_species = {}
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


def set_disamar_cfg_sim(cfg, dis_cfg, ground_points, profiles, albedo, i_t, i_x):

    # adapted from Pepijn's E2E
    # Modify the disamar input file for simulation of spectra
    # all instruments modification off, no slit, high resolution


    # GENERAL

    if cfg['rtm']['dismas_sim']:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(1)
    else:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(0)


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
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_start'].setvalue(cfg['spec_settings']['wave_start'])
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_end'].setvalue(cfg['spec_settings']['wave_end'])
    dis_cfg['INSTRUMENT','wavelength_range', 'wavelength_step'].setvalue(cfg['spec_settings']['dwave'])

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
    pt_retr.set_rawvalue(np.asarray([np.round(profiles['p'][i_t,i_x,:],3), np.round(profiles['t'][i_t,i_x,:],3), np.ones(profiles['p'][i_t,i_x,:].shape, dtype=float)]).T)

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
            print('unknown S2 albedo band')
            albedo_wvl = 440.0

    dis_cfg['SURFACE', 'wavelDependentSim', 'wavelSurfAlbedo'].setvalue(albedo_wvl)
    dis_cfg['SURFACE', 'wavelDependentSim', 'surfAlbedo'].setvalue(albedo[i_t,i_x])

    dis_cfg['SURFACE','pressure', 'surfPressureSim'].setvalue(profiles['p'][i_t,i_x, 0])
    dis_cfg['SURFACE','pressure', 'surfPressureRetr'].setvalue(profiles['p'][i_t,i_x, 0])


    return dis_cfg

def convert_atm_profiles(cfg, atm):
    # convert atmospheric profiles for usage in disamar

    profiles = {}

    if config['atm']['type'] == 'afgl':
         # afgl, input is in molec/m2, convert to ppmv

        nact, nalt, nlev =  atm[0, 0].plev.shape

        for var in cfg['atm']['gases'].copy().append('p','t'):
            profiles[var] = np.zeros((nact, nalt, nlev))

        for idx in range(nact):
            for idy in range(nalt):

                psfc = atm[i_t, i_x].psurf / 100.0 # Pa --> hPa

                plev = atm[i_t, i_x].plev[::-1] / 100.0 # levels,  Pa --> hPa, surface --> TOA
                play = atm[i_t, i_x].play[::-1] / 100.0 # layers, Pa --> hPa, surface --> TOA

                tlev = atm[i_t, i_x].tlev[::-1]         # levels,  K --> hPa, surface --> TOA
                tlay = atm[i_t, i_x].tlay[::-1]         # layers, K --> hPa, surface --> TOA

                air = atm[i_t, i_x].air[::-1]           # layers, molec/m2, surface --> TOA
                
                profiles['p'][i_t, i_x, :] = plev
                profiles['t'][i_t, i_x, :] = tlev

                for gas in cfg['atm']['gases']:
                    gas_prof = atm[i_t, i_x].NO2[::-1]           # layers, molec/m2, surface --> TOA

                    # convert molec/m2 to volume mixing ratio
                    gas_ppmv = gas_prof / air * 1e6 # molec/m2 --> ppmv

                    # interpolate layers to levels using log(p)

                    # no extrapolation, edge value is repeated
                    profiles[gas][i_t, i_x, :] = np.interp( np.log(plev[::-1]), np.log(play[::-1]), gas_ppmv[::-1] )[::-1] # strictly increasing

                    # extrapolation
                    # interp_gas = scipy.interpolate.RegularGridInterpolator(np.reshape(np.log(play.T),(1,len(play))), gas_ppmv, method='linear',bounds_error=False, fill_value=None)
                    # profiles[gas][i_t, i_x, :] = interp_gas(np.log(plev))

        # # expand profiles to TOA, otherwise disamar does not like it
        # profiles['p'] = np.append(profiles['p'],0.015)
        # profiles['t'] = np.append(profiles['t'],173.526)
        # if 'no2' in cfg['atm']['gases']:
        #     profiles['no2'] = np.append(profiles['no2'],6.2958549E-09)
        # if 'o3' in cfg['atm']['gases']:
        #     profiles['o3'] = np.append(profiles['o3'],7.4831730E-01)

    elif config['atm']['type'] == 'cams':
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
            var_layer = np.transpose(atm[var][::-1], axes=[1,2,0]) # layers, ppmv, surface --> TOA

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

def run_disamar(filename):

    dis_cfg = libRT_no2.RT_configuration(filename=filename)
    output_filename = filename.replace('.in', '.h5')

    RT = libRT_no2.rt_run(cfg=dis_cfg, disamar=config['rtm']['disamar_exe'], output=output_filename, quiet=True, debug=False)

    cwd = os.getcwd()

    try:
        print('starting: ' + filename)
        RT()
        os.chdir(cwd)
        return 0
    except:
        print('Failed: ' + filename)
        os.chdir(cwd)
        return -1

    return


# >>>>>>>>>>> to do here: move functions above to lib


def scene_generation_module_nitro(config):
    """
    Scene generation algorithm for NO2

    Note: for now only works for profile: orbit
    """

    # =============================================================================================
    # 1) get the geometry data
    # =============================================================================================

    gm_data = get_gm_data(config['input']['gm'])
    nalt, nact = gm_data['sza'].shape

    # =============================================================================================
    # 2) get collocated S2 albedo data
    # =============================================================================================

    albedo = np.zeros([nalt, nact])

    file_exists = os.path.isfile(config['S2_albedo']['dump'])
    if (file_exists and (not config['S2_albedo']['forced'])):
        albedo = np.load(config['S2_albedo']['dump'])
    else:
        print('S2 data')
        albedo = libSGM.get_sentinel2_albedo(gm_data, config,band=config['S2_albedo']['band'])
        np.save(config['S2_albedo']['dump'], albedo)    # .npy extension is added if not given
    

    # =============================================================================================
    # 3) get collocated microHH data and apply convolution
    # =============================================================================================
    if config['atm']['microHH']['use']:
        if ((not os.path.exists(config['atm']['microHH']['dump'])) or config['atm']['microHH']['forced']):
            meteodata = libATM.get_atmosphericdata(gm_data['lat'], gm_data['lon'], config['atm']['microHH'],
                                                   config['kernel_parameter'])
            # Dump microHH dictionary into temporary pkl file
            pickle.dump(meteodata.__dict__, open(config['atm']['microHH']['dump'], 'wb'))
        else:
            # Read microHH from pickle file
            meteodata = Dict2Class(pickle.load(open(config['atm']['microHH']['dump'], 'rb')))

    # =============================================================================================
    # 4) build model atmosphere and write output
    # =============================================================================================
    match config['atm']['type']:

        case 'afgl':
            # 4A) get a model atm from AFGL files

            nlay = config['atm']['afgl']['nlay']  # number of layers
            dzlay = config['atm']['afgl']['dzlay']
            # we assume the same standard atm for all pixels of the granule
            atm = libATM.get_AFGL_atm_homogenous_distribution(config['atm']['afgl']['path'], nlay, dzlay)
            
            if config['atm']['microHH']['use']:
                #combine the microHH data with standard atm afgl, regridding microHH to std atm layers
                atm = libATM.combine_meteo_standard_atm(meteodata, atm, config['atm']['microHH']['gases'])

            sgm_output_atm_afgl(config['output']['atm'], atm, albedo, gm_data, meteodata, config['atm']['gases'])

        case 'cams':
            # 4B) get atm from CAMS

            # interpolate CAMS field
            atm = get_cams_profiles(config, config['atm']['cams']['start'], gm_data)

            if config['atm']['dem']['use']:
                # correct surface pressure with DEM
                atm['zsfc'] = get_dem(config['atm']['dem']['path'], gm_data['lat'], gm_data['lon'])
                atm['psfc'] = atm['psl'] * np.exp( -1 * atm['zsfc'] / 8000.) # use 8 km scale height

            if config['atm']['microHH']['use']:
                # combine the microHH data with CAMS atm, regridding microHH to CAMS layers
                combine_mhh_cams( meteodata, atm, species=['no2'])

            sgm_output_atm_cams(config, atm, albedo, gm_data)

    # =============================================================================================
    # 5) radiative transfer simulations with DISAMAR
    # =============================================================================================

    print('radiative transfer simuation...')

    # 5A) generate the disamar config files and tmp dirs
    print('creating config files disamar')

    # convert atm profiles to disamar format
    dis_profiles = convert_atm_profiles(config, atm)

    dis_cfg = libRT_no2.RT_configuration(filename=config['rtm']['disamar_cfg_template'])
    dis_cfg_filenames=[]
    timestamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    for ialt in range(nalt):
        for iact in range(nact):

            dis_cfg = set_disamar_cfg_sim(config, dis_cfg, gm_data, dis_profiles, albedo, ialt, iact)

            tmp_dir = '{}/{}/{:03d}'.format( config['rtm']['tmp_dir'], timestamp, iact)
            if not os.path.isdir(tmp_dir):
                os.makedirs(tmp_dir, exist_ok=True)

            filename = '{}/{}_{:03d}.in'.format(tmp_dir, ialt, iact)
            dis_cfg.write(filename=filename)
            dis_cfg_filenames.append(filename)
            print(filename)

    # >>>>>>>>>> tmp
            break
        break
    # >>>>>>>>>> tmp

    # 5B) run disamar in parallel
    print('running disamar')

    with multiprocessing.Pool(config['rtm']['n_threads']) as pool:
        stat = set(pool.map(run_disamar, dis_cfg_filenames))

    # =============================================================================================
    # 6) sgm output to radiometric file
    # =============================================================================================

    # # define line-by-line wavelength grid
    # rad_output = {}
    # rad_output['wavelength_lbl'] = np.arange(config['spec_settings']['wave_start'], config['spec_settings']['wave_end'],
    #                                          config['spec_settings']['dwave'])  # nm
    # >>>>>>>>>> to do here: write output
    # sgm_output_radio(config['rad_output'], rad_output)

    print('=>sgm calculation finished successfully')
    return


if __name__ == '__main__':
    
    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )


    config = yaml.safe_load(open(sys.argv[1]))
    scene_generation_module_nitro(config)

