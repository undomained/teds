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
from scipy.interpolate import RegularGridInterpolator

from lib import libATM, libSGM, libRT_no2, libNumTools, constants
from lib.libWrite import writevariablefromname

logger = logging.getLogger('E2E')

class Emptyclass:
    """Empty class. Data container."""
    
    pass

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

def get_cams_profiles( cfg, time, lats, lons):


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
        
        if var != 't':
            var_name = 'ppmv_'+var
        elif var == 't':
            var_name = 'tlay'

        profile_data[var_name] = np.array( varlist )

        # convert mmr to ppmv
        if var in mmr_to_ppmv:
            profile_data[var_name] *= mmr_to_ppmv[var]

        # swap axes
        profile_data[var_name] = profile_data[var_name].transpose(1,2,0)

    profile_data['lat'] = lats
    profile_data['lon'] = lons


    # calculate pressure profiles
    profile_data['play'] = (1e-2 * profile_data['hyam'][np.newaxis,np.newaxis,...] + profile_data['psfc'][...,np.newaxis] * profile_data['hybm'][np.newaxis,np.newaxis,...]) # layers, hPa, TOA --> surface
    profile_data['plev'] = (1e-2 * profile_data['hyai'][np.newaxis,np.newaxis,...] + profile_data['psfc'][...,np.newaxis] * profile_data['hybi'][np.newaxis,np.newaxis,...]) # levels, hPa, TOA --> surface
    
    return profile_data


def mixingratio_to_column(config,atm):
    # ppmv to partial and total column

    # calculate partial and total column
    # pressure drop per layer [Pa]
    pdlev = (atm.plev[:,:,1:] - atm.plev[:,:,:-1])*1e2

    for gas in config['atm']['gases']:

        # ppmv profile
        gas_ppmv = atm.__getattribute__('ppmv_'+gas) # [ppmv]

        # Convert volume mixing ratio to partial column

        # partial column in [mol/m2] at layers
        gas_partialcolumn = pdlev * gas_ppmv*1e-6 / ( constants.g0 * constants.MDRYAIR ) 
        # convert to [molec/cm^2]
        gas_partialcolumn *=  constants.NA * 1e-4

        # calculate total column
        gas_totalcolumn = gas_partialcolumn.sum(axis=-1)

        atm.__setattr__('dcol_'+gas, gas_partialcolumn)
        atm.__setattr__('col_'+gas, gas_totalcolumn)
    
    return atm

def column_to_mixingratio(atm):
    # partial column to mixing ratio and total column

    atm.ppmv_no2 = atm.dcol_no2/atm.dcol_air*1.e6 # [ppmv]
    atm.ppmv_o3 = atm.dcol_o3/atm.dcol_air*1.e6 # [ppmv]

    # NO2 and O3 usually in molec/cm2
    atm.dcol_no2 *= 1e-4 # [molec/cm2]
    atm.dcol_o3 *= 1e-4 # [molec/cm2]

    atm.col_no2 = np.sum(atm.dcol_no2,axis=-1)    # [molec/cm2]
    atm.col_o3 = np.sum(atm.dcol_o3,axis=-1)   # [molec/cm2]
    
    return atm

def elevation_to_pressure(atm,psurf):
    # calculate pressures from elevation using scale height and surface pressure
    # output in [hPa]

    scale_height = 8000. #[m]
    atm.plev = psurf[...,np.newaxis]* np.exp( -1 * atm.zlev / scale_height)
    atm.play = psurf[...,np.newaxis] * np.exp( -1 * atm.zlay / scale_height)

    return atm

def combine_mhh_cams( mhh_data, atm, species=['no2', 'co2', 'no']):
    
    nalt, nact = atm.psfc.shape

    # convert mhh partial column [molec./m2] to ppmv

    # pressure drop per layer
    mhh_dplev = (mhh_data.plev[:,:,:-1] - mhh_data.plev[:,:,1:]) *1e2 # [Pa]

    for s in species:
        setattr( mhh_data, 'ppmv_'+s, getattr(mhh_data, s+'_raw') * constants.g0 * constants.MDRYAIR * 1e6 / (constants.NA * mhh_dplev.data ) ) #[ppmv]

    # take mean over mhh layers inside cams layer

    nz = atm.play.shape[-1] # output number of layers

    for s in species:
        setattr( mhh_data, 'ppmv_'+s+'_regrid', np.zeros(( nalt, nact, nz)) )

    for ix in tqdm.tqdm(range(nalt)):
        for iy in range(nact):
            for iz in range(nz):

                idx = np.where((mhh_data.play[ix,iy,:] >= atm.plev[ix,iy,iz]) & (mhh_data.play[ix,iy,:] < atm.plev[ix,iy,iz+1]))[0]

                if len(idx) < 1:
                    continue

                for s in species:
                    getattr(mhh_data, 'ppmv_'+s+'_regrid')[ix,iy,iz] = np.nanmean(getattr(mhh_data, 'ppmv_'+s)[ix,iy,idx])

    # combine mhh vertical regrid and cams
    for s in species:
        ppmv_combined = getattr(atm, 'ppmv_'+s) + getattr(mhh_data, 'ppmv_'+s+'_regrid')
        setattr(atm, 'ppmv_'+s, ppmv_combined)

    atm.__setattr__('xpos',mhh_data.x_new[:])
    atm.__setattr__('ypos',mhh_data.y_new[:])

    return atm

def add_clouds(cfg, atm):
    # Add cloud parameters to atm for specified range

    atm.cf = np.zeros_like(atm.lat)
    atm.cot = np.ma.masked_all_like(atm.lat)
    atm.cbp = np.ma.masked_all_like(atm.lat)
    atm.ctp = np.ma.masked_all_like(atm.lat)

    if 'alt' in cfg['atm']['cloud']:
        slice_alt = slice(cfg['atm']['cloud']['alt']['start'],cfg['atm']['cloud']['alt']['stop']+1)
    else:
        slice_alt = slice(0,None)
    if 'act' in cfg['atm']['cloud']:
        slice_act = slice(cfg['atm']['cloud']['act']['start'],cfg['atm']['cloud']['act']['stop']+1)
    else:
        slice_act = slice(0,None)

    atm.cf[slice_alt, slice_act] =  cfg['atm']['cloud']['cloud_fraction']
    atm.cot[slice_alt, slice_act] =  cfg['atm']['cloud']['cloud_optical_thickness']
    atm.cbp[slice_alt, slice_act] = cfg['atm']['cloud']['cloud_bottom_pressure']
    atm.ctp[slice_alt, slice_act] = cfg['atm']['cloud']['cloud_top_pressure']

    return atm


def convolvedata_nitro(atm, albedo, microhh, config):

    """Convolve meteo and albedo data.

    Parameters
    ----------
    meteodata : Class
        Meteo data
    config : Dict
        Dict containing configuration parameters.

    """
    
    dx = np.mean(microhh.dx)
    dy = np.mean(microhh.dy)
    
    conv_settings = libNumTools.getconvolutionparams(config['kernel_parameter'], dx, dy)
    
    # convolution of albedo   
    atm.__setattr__("albedo_conv", libNumTools.convolution_2d(albedo, conv_settings))
    
    for gas in config['atm']['gases']:
        dcolgas = atm.__getattribute__("dcol_"+gas)
        dcolgas_conv = np.zeros_like(dcolgas)
        ppmvgas = atm.__getattribute__("ppmv_"+gas)
        ppmvgas_conv = np.zeros_like(ppmvgas)
        for iz in range(atm.play[0,0,:].size):
            # convolution
            dcolgas_conv[:, :, iz] = libNumTools.convolution_2d(dcolgas[:, :, iz], conv_settings)
            ppmvgas_conv[:, :, iz] = libNumTools.convolution_2d(ppmvgas[:, :, iz], conv_settings)
            atm.__setattr__("dcol_"+gas+"_conv", dcolgas_conv)
            atm.__setattr__("ppmv_"+gas+"_conv", ppmvgas_conv)

    return atm


def interpolate_data_regular_nitro(indata, gm_data, config):
    """Interpolate data on regular cartesian grid.

    Meteo data is in regular x-y grid and this is used to interpolate
    the values to x-y of gm grid. This is implementation is faster
    compare to interpolating data in lat-lon grid.

    Parameters
    ----------
    indata : Class
        Meteo data
    gm_data : Dict
        Parameters from geometry module
    gases : List
        List of gases to be processed

    Returns
    -------
    albedo : Matrix
        Albedo on gm grid
    outdata: Class
        Meteo data on gm grid
    """
    logger.debug('Interpolating data to GM mesh...')
    outdata = Emptyclass()
    outdata.__setattr__("lat", gm_data['lat'])
    outdata.__setattr__("lon", gm_data['lon'])

    dim_alt, dim_act = outdata.lat.shape   # dimensions
    dim_lay = indata.play.shape[2]
    dim_lev = indata.plev.shape[2]
    
    # Interpolate albedo to GM grid
    fa = RegularGridInterpolator((indata.ypos, indata.xpos), indata.albedo_conv, 
                                 bounds_error=False, fill_value=0.0)
    albedo = fa((gm_data['ypos'], gm_data['xpos']))
    for gas in config['atm']['gases']:
        interpdata = np.zeros([dim_alt, dim_act, dim_lay])
        conv_gas = indata.__getattribute__("dcol_"+gas+"_conv")
        for iz in range(dim_lay):
            gasmin = np.min(conv_gas[:, :, iz])
            fa = RegularGridInterpolator((indata.ypos, indata.xpos), conv_gas[:, :, iz], 
                                         bounds_error=False, fill_value=gasmin)
            interpdata[:, :, iz] = fa((gm_data['ypos'], gm_data['xpos']))
        # do not allow negative values
        interpdata = interpdata.clip(min=0)
        outdata.__setattr__('dcol_'+gas, interpdata)

        interpdata = np.zeros([dim_alt, dim_act, dim_lay])
        ppmv_gas = indata.__getattribute__("ppmv_"+gas+"_conv")
        for iz in range(dim_lay):
            gasmin = np.min(ppmv_gas[:, :, iz])
            fa = RegularGridInterpolator((indata.ypos, indata.xpos), ppmv_gas[:, :, iz], 
                                         bounds_error=False, fill_value=gasmin)
            interpdata[:, :, iz] = fa((gm_data['ypos'], gm_data['xpos']))
        # do not allow negative values
        interpdata = interpdata.clip(min=0)
        outdata.__setattr__('ppmv_'+gas, interpdata)

    if config['atm']['type'] == 'afgl':
        #Here, we make a shortcut using a equally vertically gridded atmosphere 
        #So we duplicate the vertical grid!
        dum = np.tile(indata.play[0,0,:], dim_alt*dim_act)
        outdata.__setattr__("play", np.reshape(dum,(dim_alt,dim_act,dim_lay)))

        dum = np.tile(indata.plev[0,0,:], dim_alt*dim_act)
        outdata.__setattr__("plev", np.reshape(dum,(dim_alt,dim_act,dim_lev)))

        dum = np.tile(indata.tlay[0,0,:], dim_alt*dim_act)
        outdata.__setattr__("tlay", np.reshape(dum,(dim_alt,dim_act,dim_lay)))

        dum = np.tile(indata.tlev[0,0,:], dim_alt*dim_act)
        outdata.__setattr__("tlev", np.reshape(dum,(dim_alt,dim_act,dim_lev)))

    elif config['atm']['type'] == 'cams':
        varlist = ['plev','play','tlay']
        for var in varlist:
            if 'lev' in var:
                dim_z = dim_lev
            elif 'lay' in var:
                dim_z = dim_lay
            interpdata = np.zeros([dim_alt, dim_act, dim_z])
            vardata = indata.__getattribute__(var)
            for iz in range(dim_z):
                varmin = np.min(vardata[:, :, iz])
                fa = RegularGridInterpolator((indata.ypos, indata.xpos), vardata[:, :, iz], 
                                            bounds_error=False, fill_value=varmin)
                interpdata[:, :, iz] = fa((gm_data['ypos'], gm_data['xpos']))
            # do not allow negative values
            interpdata = interpdata.clip(min=0)
            outdata.__setattr__(var, interpdata)
     
    return albedo, outdata

def recalc_total_column(atm, config):
    # recalculate total column from partial columns after convolution and regridding

    for gas in config['atm']['gases']:
        dcol = atm.__getattribute__("dcol_"+gas)

        col = np.sum(dcol, axis=-1)

        atm.__setattr__('col_'+gas, col)

    return atm


def convert_atm_to_disamar(atm, cfg):
    # convert atmospheric profiles and parameters for usage in disamar

    extrapolate_to_surface = True

    atm_disamar = {}
        
    nalt, nact, nlev =  atm.plev.shape

    variables = cfg['atm']['gases'].copy()
    variables.extend(['t','p'])
    for var in variables:
        atm_disamar[var] = np.zeros((nalt, nact, nlev))

    atm_disamar['p'] = atm.plev[:,:,::-1] # levels, hPa, surface --> TOA

    # gases
    for gas in cfg['atm']['gases']:
        gas_ppmv = atm.__getattribute__('ppmv_'+gas)   # layers, ppmv, TOA --> surface

        for idx in range(nalt):
            for idy in range(nact):

                # interpolate layers to levels using log(p)

                if extrapolate_to_surface:
                    # extrapolation
                    interp_gas = RegularGridInterpolator(np.reshape(np.log(atm.play[idx,idy,:].T),(1,len(atm.play[idx,idy,:]))), gas_ppmv[idx,idy,:], method='linear',bounds_error=False, fill_value=None)
                    atm_disamar[gas][idx, idy, :] = interp_gas(np.log(atm.plev[idx,idy,:]))[::-1]
                else:
                    # no extrapolation, edge value is repeated
                    atm_disamar[gas][idx, idy, :] = np.interp( np.log(atm.plev[idx,idy,:]), np.log(atm.play[idx,idy,:]), gas_ppmv[idx,idy,:] )[::-1] # strictly increasing
    
    # temperature
    if cfg['atm']['type'] == 'afgl':
        atm_disamar['t'] = atm.tlev[:,:,::-1] # levels,  K , surface --> TOA

    elif cfg['atm']['type'] == 'cams':
        tlay = atm.__getattribute__('tlay')   # layers, K, TOA --> TOA
        for idx in range(nalt):
            for idy in range(nact):
                if extrapolate_to_surface:
                    # extrapolation
                    interp_gas = RegularGridInterpolator(np.reshape(np.log(atm.play[idx,idy,:].T),(1,len(atm.play[idx,idy,:]))), tlay[idx,idy,:], method='linear',bounds_error=False, fill_value=None)
                    atm_disamar['t'][idx, idy, :] = interp_gas(np.log(atm.plev[idx,idy,:]))[::-1]
                else:
                    # no extrapolation, edge value is repeated
                    atm_disamar['t'][idx, idy, :] = np.interp( np.log(atm.plev[idx,idy,:]), np.log(atm.play[idx,idy,:]), tlay[idx,idy,:] )[::-1] # strictly increasing

        # omit TOA layer CAMS, pressure is 0, disamar does not like it
        atm_disamar['p'] = atm_disamar['p'][:,:,:-1]
        atm_disamar['t'] = atm_disamar['t'][:,:,:-1]
        for gas in cfg['atm']['gases']:
            atm_disamar[gas] = atm_disamar[gas][:,:,:-1]

    # set last layer of profiles to TOA, otherwise disamar does not like it
    if (atm_disamar['p'][:,:,-1] > 0.3 ).any():
        atm_disamar['p'][:,:,-1] = 0.3
        atm_disamar['t'][:,:,-1] = 250.0
        if 'no2' in cfg['atm']['gases']:
            atm_disamar['no2'][:,:,-1] = 6.2958549E-09
        if 'o3' in cfg['atm']['gases']:
            atm_disamar['o3'][:,:,-1] = 7.4831730E-01

    # do not allow negative or zero values
    for key in atm_disamar:
        atm_disamar[key] = atm_disamar[key].clip(min=1.0e-10)

    # O2-O2
    # note that for a collision complex the parent gas has to be specified here, 
    # e.g. for O2-O2 the volume mixing ratio of O2 has to be specified

    # mixing ratio of O2 is 0.20946 taken from R. Goody, Principles of atmospheric physics and chemistry,
    # Table 1.2, Oxford University Press, New York, 1995. [DAK uses 0.209476  (US Stand. Atm., 1976)]
    o2_mixing_ratio =  20.94600E+04
    if cfg['rtm']['o2o2']:
        cfg['atm']['gases'].append('o2-o2')
        atm_disamar['o2-o2'] = np.ones_like(atm_disamar['p'])*o2_mixing_ratio

    if cfg['atm']['cloud']['use']:
        atm_disamar['cloud_fraction'] = atm.cf.copy()
        atm_disamar['cloud_optical_thickness'] = atm.cot.copy()
        atm_disamar['cloud_top_pressure'] = atm.ctp.copy()
        atm_disamar['cloud_bottom_pressure'] = atm.cbp.copy()

    return atm_disamar



def set_disamar_cfg_sim(cfg, dis_cfg, ground_points, atm_disamar, albedo, i_t, i_x):

    # adapted from Pepijn's E2E
    # Modify the disamar input file for simulation of spectra
    # all instruments modification off, no slit, high resolution


    # GENERAL

    if cfg['rtm']['dismas_sim']:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(1)
        dis_cfg['GENERAL','method', 'ignoreSlitSim'].setvalue(0)
        dis_cfg['INSTRUMENT','slit_parameters', 'FWHM_irradiance_sim'].setvalue(cfg['rtm']['dismas_fwhm'])
        dis_cfg['INSTRUMENT','slit_parameters', 'FWHM_radiance_sim'].setvalue(cfg['rtm']['dismas_fwhm'])

    else:
        dis_cfg['GENERAL','method', 'simulationMethod'].setvalue(0)
        dis_cfg['GENERAL','method', 'ignoreSlitSim'].setvalue(1)

    dis_cfg['GENERAL','overall', 'numberTraceGases'].setvalue(len(cfg['atm']['gases']))

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
    pt_sim.set_rawvalue(np.asarray([atm_disamar['p'][i_t,i_x,:], atm_disamar['t'][i_t,i_x,:]]).T)
    pt_retr = dis_cfg['PRESSURE_TEMPERATURE', 'PT_retr', 'PT']
    pt_retr.set_rawvalue(np.asarray([atm_disamar['p'][i_t,i_x,:], atm_disamar['t'][i_t,i_x,:], np.ones(atm_disamar['p'][i_t,i_x,:].shape, dtype=float)]).T)

    # gas profiles
    for gas in cfg['atm']['gases']:
        gas_vmr = np.asarray([atm_disamar['p'][i_t,i_x,:], atm_disamar[gas][i_t,i_x,:]]).T
        gas_vmr_error = np.asarray([atm_disamar['p'][i_t,i_x,:], atm_disamar[gas][i_t,i_x,:], np.ones(atm_disamar['p'][i_t,i_x,:].shape, dtype=float) * 20.]).T
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

    dis_cfg['SURFACE','pressure', 'surfPressureSim'].setvalue(atm_disamar['p'][i_t,i_x, 0])
    dis_cfg['SURFACE','pressure', 'surfPressureRetr'].setvalue(atm_disamar['p'][i_t,i_x, 0])


    # clouds (HG scattering)
    if atm_disamar['cloud_fraction'][i_t,i_x] > 0.0:

        # cloud fraction
        dis_cfg['CLOUD_AEROSOL_FRACTION','wavelIndependentSim', 'fraction'].setvalue( atm_disamar['cloud_fraction'][i_t,i_x] )
        
        # cloud optical thickness
        dis_cfg['CLOUD','HGscatteringSim', 'opticalThickness'].setvalue( [ 2.0 , atm_disamar['cloud_optical_thickness'][i_t,i_x] ] )

        # cloud bottom and top pressures
        dis_cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures', 'topPressureSim'].setvalue( [ atm_disamar['cloud_bottom_pressure'][i_t,i_x], atm_disamar['cloud_top_pressure'][i_t,i_x], 0.30] )

        # radiative transfer settings for optically thick cloud
        dis_cfg['RADIATIVE_TRANSFER','numDivPointsAlt', 'numDivPointsAltSim'].setvalue( [8, int(np.ceil(1.5*atm_disamar['cloud_optical_thickness'][i_t,i_x])), 32] )
        dis_cfg['RADIATIVE_TRANSFER','RTM_Sim_Retr', 'useAddingSim'].setvalue(1)
        dis_cfg['RADIATIVE_TRANSFER','RTM_Sim_Retr', 'nstreamsSim'].setvalue(64)


    return dis_cfg



def run_disamar(filename,disamar_exe):


    dis_cfg = libRT_no2.RT_configuration(filename=filename)
    output_filename = filename.replace('.in', '.h5')

    if logger.isEnabledFor(logging.DEBUG):
        quiet = False
        debug = True
    else:
        quiet = True
        debug = False

    try:
        RT = libRT_no2.rt_run(cfg=dis_cfg, disamar=disamar_exe, output=output_filename, quiet=quiet, debug=debug)
        starttime = time.time()
        RT()
        logger.debug(f'finished: {filename} in {np.round(time.time()-starttime,1)} s')
        return 0
    except:
        logger.error(f'failed: {filename}')
        return -1


def read_disamar_output(gm_data,tmp_dir):

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
    output_rad = nc.Dataset(config['io']['sgm_rad'], mode='w')
    output_rad.title = 'Tango E2ES SGM radiometric scene'
    # output_rad.config = str(config)
    output_rad.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    output_rad.product_type = 'SGM'

    dims = ('along_track', 'across_track', 'wavelength')

    output_rad.createDimension(dims[0], nalt)     # along track axis
    output_rad.createDimension(dims[1], nact)     # across track axis
    output_rad.createDimension(dims[2], nlbl)     # spectral axis

    # grp = output_rad.createGroup('science_data')

    # wavelength
    _ = writevariablefromname(output_rad, 'wavelength', ('wavelength',), rad_output['wavelength_lbl'])
    # solar irradiance
    _ = writevariablefromname(output_rad, 'solarirradiance', ('wavelength',), rad_output['solar_irradiance'])
    
    # radiance
    _ = writevariablefromname(output_rad, 'radiance_sgm', dims, rad_output['radiance'])
    
    output_rad.close()



def sgm_output_atm(config, atm, albedo, microhh_data, mode='raw'):
    # write atm to netcdf
    # mode = {'raw' , 'convolved'}

    gases = config['atm']['gases']

    nalt, nact = atm.lat.shape
    nlay, nlev = atm.play.shape[-1], atm.plev.shape[-1]

    # create file
    if mode == 'raw':
        comment = 'Raw, non-convolved scene'
        file = config['io']['sgm_atm_raw']
    elif mode == 'convolved':
        comment = 'Instrument point-spead function convolved scene'
        file = config['io']['sgm_atm']

    output_atm = nc.Dataset(file, mode='w')
    output_atm.title = 'Tango E2ES Nitro SGM atmospheric scene'
    # output_atm.config = str(config)
    output_atm.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    output_atm.comment = comment



    output_atm.createDimension('along_track', nalt)      # along track axis
    output_atm.createDimension('across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis
    output_atm.createDimension('location', 3)                 # source coordiantes
    output_atm.createDimension('emission', 1)                 # emission strength

    # add variables
    dims_layer = ('along_track', 'across_track', 'number_layers')
    dims_level = ('along_track', 'across_track', 'number_levels')
    dims_2d = ('along_track', 'across_track')

    _ = writevariablefromname(output_atm, 'latitude', dims_2d, atm.lat)
    _ = writevariablefromname(output_atm, 'longitude', dims_2d, atm.lon)
    _ = writevariablefromname(output_atm, 'pressure_levels', dims_level, atm.plev)
    _ = writevariablefromname(output_atm, 'pressure_layers', dims_layer, atm.play)

    if hasattr(atm, 'tlay'):
        _ = writevariablefromname(output_atm, 'temperature', dims_layer, atm.tlay)

    if mode == 'raw':
        bands = [x.removesuffix('_albedo') for x in albedo.__dict__.keys() if '_albedo' in x]
        for band in bands:
            var_alb = writevariablefromname(output_atm, 'albedo_'+ band, dims_2d, 
                                            albedo.__getattribute__(band+'_albedo'))
            var_alb.setncattr("central wavelength", albedo.__getattribute__(band+'_central_wave'))
            var_alb.setncattr("band width", albedo.__getattribute__(band+'_band_width'))
    elif mode == 'convolved':
        var_alb = writevariablefromname(output_atm, 'albedo', dims_2d, albedo)

    for gas in gases:
        _ = writevariablefromname(output_atm, 'conc_'+ gas, dims_layer, atm.__getattribute__('ppmv_'+gas)) # [ppmv]
        _ = writevariablefromname(output_atm, 'subcol_density_'+gas, dims_layer, atm.__getattribute__('dcol_'+gas)) # [molec/cm2]
        _ = writevariablefromname(output_atm, 'column_'+gas, dims_2d, atm.__getattribute__('col_'+gas)) # [molec/cm2]

    # clouds
    if config['atm']['cloud']['use'] and mode=='convolved':
        _ = writevariablefromname(output_atm, 'cloud_fraction', dims_2d, atm.cf)
        _ = writevariablefromname(output_atm, 'cloud_optical_thickness', dims_2d, atm.cot)
        _ = writevariablefromname(output_atm, 'cloud_bottom_pressure', dims_2d, atm.cbp)
        _ = writevariablefromname(output_atm, 'cloud_top_pressure', dims_2d, atm.ctp)



    if config['atm']['microHH']['use'] and mode=='raw':
        #information on emission source    
        substr = 'source'
        attr = microhh_data.__dict__.keys()
        attr_src = [string for string in attr if substr in string]
        for src in attr_src:
            _ = writevariablefromname(output_atm, src+'_location', 'location', microhh_data.__getattribute__(src))

        substr = 'emission'
        attr = microhh_data.__dict__.keys()
        attr_emi = [string for string in attr if substr in string]
        for emi in attr_emi:
            _ = writevariablefromname(output_atm, emi.removesuffix('_in_kgps'), 'emission' , microhh_data.__getattribute__(emi))

        # xpos and ypos       
        _ = writevariablefromname(output_atm, 'xpos', 'across_track', microhh_data.x_new)
        _ = writevariablefromname(output_atm, 'ypos', 'along_track', microhh_data.y_new)
        # # level height
        # _ = writevariablefromname(output_atm, 'levelheight', dims_level, atm.zlev)
        # # central layer height
        # _ = writevariablefromname(output_atm, 'central_layer_height', dims_layer, atm.zlay)

    output_atm.close()

    return


def scene_generation_module_nitro(config):
    """
    Scene generation algorithm for NO2

    Note: for now only works for profile: orbit
    """

    logger.info('Starting SGM calculation')

    start_time = time.time()

    # =============================================================================================
    # 1) get the geometry data
    # =============================================================================================

    gm_data = get_gm_data(config['io']['gm'])
    nalt, nact = gm_data['sza'].shape
    logger.info(f"Pixels along track: {nalt} ; Pixels across track: {nact}")

    # =============================================================================================
    # 2) get microHH data
    # =============================================================================================
    if config['atm']['microHH']['use']:
        if ((not os.path.exists(config['atm']['microHH']['dump'])) or config['atm']['microHH']['forced']):
            logger.info(f"Loading microHH data: {config['atm']['microHH']['path_data']}")

            microhh_data = libATM.get_atmosphericdata_new(gm_data['lat'], gm_data['lon'], config['atm']['microHH'])
            
            # Dump microHH dictionary into temporary pkl file
            pickle.dump(microhh_data.__dict__, open(config['atm']['microHH']['dump'], 'wb'))
        else:
            logger.info(f"Loading microHH data from dump file: {config['atm']['microHH']['dump']}")
            # Read microHH from pickle file
            microhh_data = Dict2Class(pickle.load(open(config['atm']['microHH']['dump'], 'rb')))
        lat = microhh_data.lat
        lon = microhh_data.lon
    else:
        microhh_data = Emptyclass()
        lat = gm_data['lat']
        lon = gm_data['lon']

    # =============================================================================================
    # 3) get collocated S2 albedo data on microHH grid
    # =============================================================================================

    albedo = np.zeros([nalt, nact])

    file_exists = os.path.isfile(config['S2_albedo']['dump'])
    if (file_exists and (not config['S2_albedo']['forced'])):
        logger.info(f"Loading S2 data from dump file: {config['S2_albedo']['dump']}")
        albedo = pickle.load(open(config['S2_albedo']['dump'], 'rb'))
    else:
        logger.info(f"Downloading S2 data for band: {config['S2_albedo']['band']}")
                
        albedo = libSGM.get_sentinel2_albedo_new(lat, lon, [config['S2_albedo']['band']])

        # clip albedo
        bands = [x for x in albedo.__dict__.keys() if '_albedo' in x]
        for band in bands:
            albedo_tmp = albedo.__getattribute__(band)
            albedo_tmp[albedo_tmp>1.0] = 1.0
            albedo_tmp[albedo_tmp<0.0] = 0.0
            albedo.__setattr__(band,albedo_tmp)
        
        pickle.dump(albedo, open(config['S2_albedo']['dump'], 'wb')) # .pkl extension is added if not given
    
    # =============================================================================================
    # 4) build model atmosphere and write output
    # =============================================================================================
    match config['atm']['type']:

        case 'afgl':

            logger.info(f"Loading afgl atmosphere: {config['atm']['afgl']['path']}")

            # 4A) get a model atm from AFGL files

            nlay = config['atm']['afgl']['nlay']  # number of layers
            dzlay = config['atm']['afgl']['dzlay']
            # we assume the same standard atm for all pixels of the granule
            atm = libATM.get_AFGL_atm_homogenous_distribution(config['atm']['afgl']['path'], nlay, dzlay)
            
            config['atm']['afgl'] =  {}
            config['atm']['meteo'] = {}
            config['atm']['meteo']['gases'] = config['atm']['microHH']['gases']
            config['atm']['afgl']['gases'] = config['atm']['gases']
            microhh_data.__setattr__("albedo", albedo)

            if config['atm']['microHH']['use']:
                config['atm']['afgl']['afgl_only'] =  False
            else:
                config['atm']['afgl']['afgl_only'] =  True
                microhh_data.__setattr__("lat", lat)
                microhh_data.__setattr__("lon", lon)
                microhh_data.__setattr__("x_new", np.zeros_like(lat))
                microhh_data.__setattr__("y_new", np.zeros_like(lat))

            atm = libATM.combine_meteo_standard_atm_new(microhh_data, atm, config['atm'])

            atm = column_to_mixingratio(atm)
            atm = elevation_to_pressure(atm,atm.psurf*1e-2)
            
        case 'cams':

            logger.info(f"Loading CAMS atmosphere: {config['atm']['cams']['path']}")

            # 4B) get atm from CAMS

            # interpolate CAMS field
            atm = get_cams_profiles(config, config['atm']['cams']['start'], lat, lon)
            atm = Dict2Class(atm)
            if config['atm']['dem']['use']:

                logger.info(f"Loading DEM: {config['atm']['dem']['path']}")
                # correct surface pressure with DEM
                atm.zsfc = get_dem(config['atm']['dem']['path'], lat, lon)
                atm.psfc = atm.psl * np.exp( -1 * atm.zsfc / 8000.) # use 8 km scale height

            if config['atm']['microHH']['use']:

                logger.info(f"Merging CAMS and microHH atmosphere")

                microhh_data = elevation_to_pressure(microhh_data,atm.psfc)

                # combine the microHH data with CAMS atm, regridding microHH to CAMS layers
                atm = combine_mhh_cams( microhh_data, atm, species=['no2'])
            
            atm = mixingratio_to_column(config,atm)


    logger.info(f"Writing raw scene atmosphere: {config['io']['sgm_atm_raw']}")
        
    sgm_output_atm(config, atm, albedo, microhh_data, mode = 'raw')

    # =============================================================================================
    # 5) convolve and regrid atmosphere to instrument grid
    # =============================================================================================

    # only convolve + interpolate when microHH is used
    # otherwise source grid is already interpolated to instrument grid
    if config['atm']['microHH']['use']:
        logger.info('Convolving SGM atmosphere with instrument spatial response function')

        #convolution of albeod and microHH data with instrument spatial response
        atm = convolvedata_nitro(atm, albedo.B01_albedo, microhh_data, config)
        
        # create a transform method 
        trans = libNumTools.TransformCoords(microhh_data.no2_source[1:])

        # convert lat-lon of gm to x-y and get bounds
        gm_data_class = Dict2Class(gm_data)
        gm_data['xpos'], gm_data['ypos'] = trans.latlon2xymts(gm_data_class.lat, gm_data_class.lon)

        #interpolate convolved data to gm grid
        albedo, atm = interpolate_data_regular_nitro(atm, gm_data, config)

        # recalculate total columns
        atm = recalc_total_column(atm, config)

        # add clouds
        if config['atm']['cloud']['use']:
            atm = add_clouds(config,atm)

        # write atm to file
        logger.info(f"Writing convolved scene atmosphere: {config['io']['sgm_atm']}")
        sgm_output_atm(config, atm, albedo, microhh_data, mode = 'convolved')

        
    # =============================================================================================
    # 6) radiative transfer simulations with DISAMAR
    # =============================================================================================

    logger.info('Radiative transfer simulation')

    # generate the disamar config files and tmp dirs
    logger.info('Creating config files DISAMAR')

    # convert atm profiles to disamar format
    atm_disamar = convert_atm_to_disamar(atm, config)


    timestamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    tmp_dir = '{}/{}'.format( config['rtm']['tmp_dir'], timestamp)
    if not os.path.isdir(tmp_dir):
        logger.info(f'Creating tmp directory DISAMAR: {tmp_dir}')

        os.makedirs(tmp_dir, exist_ok=True)


    if os.path.isfile(config['rtm']['disamar_cfg_template']):
        tmp_template_disamar = libRT_no2.disamar_add_gas_cfg(config, tmp_dir)
    else:
        logger.error(f"File {config['rtm']['disamar_cfg_template']} not found")

    dis_cfg_filenames=[]



    for ialt in range(nalt):
        for iact in range(nact):
            dis_cfg_template = libRT_no2.RT_configuration(filename=tmp_template_disamar)
            dis_cfg = set_disamar_cfg_sim(config, dis_cfg_template, gm_data, atm_disamar, albedo, ialt, iact)

            filename = f'{tmp_dir}/alt{ialt:04d}_act{iact:03d}.in'
            dis_cfg.write(filename=filename)
            dis_cfg_filenames.append(filename)
    

    # run disamar in parallel
    logger.info('Running DISAMAR')

    with multiprocessing.Pool(config['rtm']['n_threads']) as pool:
        # stat = pool.starmap(run_disamar, zip(dis_cfg_filenames, repeat(config['rtm']['disamar_exe'])),chunksize=1)
        stat = pool.starmap(run_disamar, tqdm.tqdm(zip(dis_cfg_filenames, repeat(config['rtm']['disamar_exe'])),total=len(dis_cfg_filenames)),chunksize=1)

    if sum(stat) != 0:
        logger.error('Error in at least one RT calculation run with DISAMAR')

    # read disamar output
    logger.info('Reading DISAMAR output')

    dis_output = read_disamar_output(gm_data, tmp_dir)

    # cleanup
    if config['rtm']['cleanup']:
        logger.info('Cleaning up tmp directory DISAMAR')
        shutil.rmtree(tmp_dir)

    # =============================================================================================
    # 7) sgm output to radiometric file
    # =============================================================================================
    logger.info(f"Writing scene radiance to: {config['io']['sgm_rad']}")

    sgm_output_radio(config, dis_output)

    logger.info(f'SGM calculation finished in {np.round(time.time()-start_time,1)} s')
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
    scene_generation_module_nitro(config)

