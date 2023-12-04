#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sept 27 17:30:22 2023
simplified level2 module
@author: manugv
"""
import sys
from netCDF4 import Dataset
from numpy import random, ones_like
import numpy as np
from ..lib.libWrite import writevariablefromname
from .l1bl2 import write_gasdata, write_proxygasdata

class DataCont:
    pass

def readsgmatmosphere(sgmatmosphere_file, instrument):
    """Read the SGM data.

    Read actual integrated gas data, lat, lon, zlay and aircolumn.

    Parameters
    ----------
    sgmatmosphere_file : String
        SGM filename.
    instrument : String
        TANGO_Carbon or TANGO_Nitro
    """
    
    if(instrument not in ["TANGO_Carbon", "TANGO_Nitro"]):
        sys.exit('wrong instrument specification in simplified L2, choose between "TANGO_Carbon" and "TANGO_Nitro"')

    data = DataCont()
    
    # read the data
    f = Dataset(sgmatmosphere_file, "r")

    data.__setattr__("zlay", f["zlay"][:])
    data.__setattr__("lon", f["lon"][:])
    data.__setattr__("lat", f["lat"][:])

    if(instrument == 'TANGO_Carbon'):
        data.__setattr__("CO2_column", f["XCO2"][:])
        data.__setattr__("CH4_column", f["XCH4"][:])
        data.__setattr__("H2O_column", f["XH2O"][:])
        data.__setattr__("CO2_subcolumn", f["dcol_co2"][:])
        data.__setattr__("CH4_subcolumn", f["dcol_ch4"][:])
        data.__setattr__("H2O_subcolumn", f["dcol_h2o"][:])
        data.__setattr__("air_column", f["col_air"][:])

    if(instrument == 'TANGO_Nitro'):
        data.__setattr__("NO2_column", f["col_co2"][:])
        data.__setattr__("NO2_subcolumn", f["dcol_co2"][:])
        data.__setattr__("air_column", f["col_air"][:])

    return data

def getlevel2retrieval(conc, precision, seed):
    """Generate noisy data.

    Parameters
    ----------
    conc : Array (m,n)
        Concentration
    precision : Array (m, n)
        Level 2 precision
    seed : Integer
        Seed for random generation

    """
    noise = random.default_rng(seed).normal(0, 1, conc.shape)
    return conc + noise*precision

def writedata_l2data(level2_outputfile, data, instrument):
    """Write level 2 data.

    Parameters
    ----------
    level2_outputfile : String
        Output file name of the Level 2 data.
    data : Data Class
        Class containing data of SGM and the newly generated level 2 data.
    gas : String
        Kind of gas. Should be in lower-case.

    """
    # filename
    output_l2 = Dataset(level2_outputfile, mode='w')
    
    if(instrument=='TANGO_Carbon'):
        output_l2.title = 'Tango Carbon E2ES L2 synthetic product'

        # create dimensions
        nalt, nact, nlay = data.zlay.shape
        output_l2.createDimension('bins_along_track', nalt)      # along track axis
        output_l2.createDimension('bins_across_track', nact)     # across track axis
        output_l2.createDimension('number_layers', nlay)         # layer axis

        # dimensions of two types
        _dims2d = ('bins_along_track', 'bins_across_track')
        _dims3d = ('bins_along_track', 'bins_across_track', 'number_layers',)

        # convergence
        _ = writevariablefromname(output_l2, 'convergence', _dims2d,  data.conv)
        # latitude
        _ = writevariablefromname(output_l2, 'latitude', _dims2d, data.lat)
        # longitude
        _ = writevariablefromname(output_l2, 'longitude', _dims2d, data.lon)
        # maxiterations
        _ = writevariablefromname(output_l2, 'maxiterations', _dims2d, data.max_iter)
        # chi2
        _ = writevariablefromname(output_l2, 'spectralchi2', _dims2d, data.chi2)
        # albedo
        _ = writevariablefromname(output_l2, 'albedo', _dims2d, data.albedo)
        # XCO2
        _ = writevariablefromname(output_l2, 'XCO2', _dims2d, data.XCO2)
        _ = writevariablefromname(output_l2, 'avgkernelXCO2', _dims3d, data.col_avg_kernel_XCO2)
        _ = writevariablefromname(output_l2, 'precisionXCO2', _dims2d, data.precision_XCO2)
        # XCH4
        _ = writevariablefromname(output_l2, 'XCH4', _dims2d, data.XCH4)
        _ = writevariablefromname(output_l2, 'avgkernelXCH4', _dims3d, data.col_avg_kernel_XCH4)
        _ = writevariablefromname(output_l2, 'precisionXCH4', _dims2d, data.precision_XCH4)
        # XH2O
        _ = writevariablefromname(output_l2, 'XH2O', _dims2d, data.XCH4)
        _ = writevariablefromname(output_l2, 'avgkernelXH2O', _dims3d, data.col_avg_kernel_XH2O)
        _ = writevariablefromname(output_l2, 'precisionXH2O', _dims2d, data.precision_XH2O)
        # XCO2 proxy
        _ = writevariablefromname(output_l2, 'proxyXCO2', _dims2d, data.XCO2_proxy)
        _ = writevariablefromname(output_l2, 'precisionproxyXCO2', _dims2d, data.precision_XCO2_proxy)
        # XCH4 proxy
        _ = writevariablefromname(output_l2, 'proxyXCH4', _dims2d, data.XCH4_proxy)
        _ = writevariablefromname(output_l2, 'precisionproxyXCH4', _dims2d, data.precision_XCH4_proxy)
        # 
        # central layer height
        _ = writevariablefromname(output_l2, 'central_layer_height', _dims3d, data.zlay)

        output_l2.close()

    if(instrument=='TANGO_Nitro'):
        output_l2.title = 'Tango Nitro E2ES L2 synthetic product'

        # create dimensions
        nalt, nact = data.lat.shape()
        nlay = data.zlay.size
        output_l2.createDimension('bins_along_track', nalt)      # along track axis
        output_l2.createDimension('bins_across_track', nact)     # across track axis
        output_l2.createDimension('number_layers', nlay)         # layer axis

        # dimensions of two types
        _dims2d = ('bins_along_track', 'bins_across_track')
        _dims3d = ('bins_along_track', 'bins_across_track', 'number_layers',)

        # convergence
        _ = writevariablefromname(output_l2, 'convergence', _dims2d,  data.conv)
        # latitude
        _ = writevariablefromname(output_l2, 'latitude', _dims2d, data.lat)
        # longitude
        _ = writevariablefromname(output_l2, 'longitude', _dims2d, data.lon)
        # maxiterations
        _ = writevariablefromname(output_l2, 'maxiterations', _dims2d, data.max_iter)
        # chi2
        _ = writevariablefromname(output_l2, 'spectralchi2', _dims2d, data.chi2)
        # albedo
        _ = writevariablefromname(output_l2, 'albedo', _dims2d, data.albedo)
        # NO2
        _ = writevariablefromname(output_l2, 'XCO2', _dims2d, data.NO2)
        _ = writevariablefromname(output_l2, 'avgkernelNO2', _dims3d, data.col_avg_kernel_NO2)
        _ = writevariablefromname(output_l2, 'precisionNO2', _dims2d, data.precision_NO2)

        output_l2.close()

def simplified_level2(config):
        
#        sgmatmosphere_inputfile, level2_outputfile, instrument, precision_rel, precision_const, seed=10):
    """Generate synthetic Level 2 data.

    Read SGM file (actual data) and create a synthetic level 2 data
    based on precision.

    Parameters
    ----------
    config['sgm_input'] : String
        SGM input file name.
    config['sl2_output'] : String
        Output file name of the Level 2 data.
    config['instrument'] : String Instrument TANGO-Carbon or TANGO-Nitro
    config['precision_relative'] : Float, relative error
    config['precision_constant'] : Float, constant error
    config['seed'] : Int, Optional
        A seed to generate random numbers for level 2 precision.

    """
    # read data for a given gas
    data = readsgmatmosphere(config['sgm_input'], config['instrument'])
    nalt, nact, nlay = data.zlay.shape

    data.__setattr__("albedo",  np.ma.MaskedArray(np.zeros((nalt,nact),np.float64), mask=True))
    data.__setattr__("chi2", np.ma.MaskedArray(np.zeros((nalt,nact), np.float64), mask=True) )
    data.__setattr__("conv", np.ma.MaskedArray(np.zeros((nalt,nact), np.int32), mask=True) )
    data.__setattr__("max_iter", np.ma.MaskedArray(np.zeros((nalt,nact), np.int32), mask=True))

    if(config['instrument']=='TANGO_Carbon'):
        data.__setattr__("precision_XCO2_proxy", data.CO2_column*config['precision relative'] + config['precision constant'])
        data.__setattr__("precision_XCH4_proxy", data.CH4_column*config['precision relative'] + config['precision constant'])
        data.__setattr__("XCO2_proxy", getlevel2retrieval(data.CO2_column, data.precision_XCO2_proxy, config['seed']))
        data.__setattr__("XCH4_proxy", getlevel2retrieval(data.CH4_column, data.precision_XCH4_proxy, config['seed']))
        data.__setattr__("col_avg_kernel_XCO2", np.ones((nalt,nact,nlay)))  # is set to one
        data.__setattr__("col_avg_kernel_XH2O", np.ones((nalt,nact,nlay)))  # is set to one
        data.__setattr__("col_avg_kernel_XCH4", np.ones((nalt,nact,nlay)))  # is set to one
        data.__setattr__("precision_XCO2", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
        data.__setattr__("precision_XCH4", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
        data.__setattr__("precision_XH2O", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
        data.__setattr__("XCO2", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
        data.__setattr__("XCH4", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
        data.__setattr__("XH2O", np.ma.MaskedArray(np.zeros((nalt,nact),), mask=True))
              
    # save the level 2 data
    writedata_l2data(config['l2_output'], data, config['instrument'])
