#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sept 27 17:30:22 2023

@author: manugv
"""

from netCDF4 import Dataset
from numpy import random, ones_like
from ..lib.libWrite import writevariablefromname


class DataCont:
    pass


def readsgmatmosphere(sgmatmosphere_file, gas):
    """Read the SGM data.

    Read actual integrated gas data, lat, lon, zlay and aircolumn.

    Parameters
    ----------
    sgmatmosphere_file : String
        SGM filename.
    gas : String
        Gas names in small caps
    """
    # create a container
    data = DataCont()
    # read the data
    f = Dataset(sgmatmosphere_file, "r")
    # read gas data
    if gas == "co2":
        gascol = f["col_co2"][:].data
    elif gas == "ch4":
        gascol = f["col_ch4"][:].data
    elif gas == "no2":
        gascol = f["col_no2"][:].data
    # read air data
    air = f["col_air"][:].data
    # read grid data
    for ky in ["lat", "lon", "zlay"]:
        if ky in f.variables.keys():
            data.__setattr__(ky, f[ky][:])
    f.close()
    data.__setattr__("actual_column", gascol)
    data.__setattr__("actual_column_air", air)
    return data


def getlevel2retrieval(conc, precision, seed):
    """Generate random data.

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


def writegasandproxydata(gas, output_l2, _dims, _dims3d, level2data, level2precision, avgkernel):
    """Write level 2 data gas and proxy data.

    Parameters
    ----------
    gas : String
        Kind of gas. Should be in lower-case.
    output_l2 : NetCDF group
        Group definition in netCDF4 file
    _dims : Tuple (String)
        Dimensions names in 2d (for netcdf)
    _dims3d : tuple (String)
        Dimensions names in 3d (for netcdf)
    level2data : Array
        Generated level2 data.
    level2precision : Array
        Generated level2 precision.
    avgkernel : Array
        Generated averaging kernel. Defaults to 1.

    """
    prec_varname = 'precision' + gas
    avgker_varname = 'avgkernel' + gas
    _ = writevariablefromname(output_l2, gas, _dims, level2data)
    _ = writevariablefromname(output_l2, prec_varname, _dims, level2precision)
    _ = writevariablefromname(output_l2, avgker_varname, _dims3d, avgkernel)
    # write data
    vargas = "proxy"+gas
    varprecgas = "precision"+vargas
    _ = writevariablefromname(output_l2, vargas, _dims, level2data)
    _ = writevariablefromname(output_l2, varprecgas, _dims, level2precision)


def writedata_l2data(level2_outputfile, data, gas):
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
    output_l2.title = 'Tango Carbon E2ES L2 synthetic product'

    # create dimensions
    nalt, nact = data.lat.shape()
    nlay = data.zlay.size
    output_l2.createDimension('bins_along_track', nalt)      # along track axis
    output_l2.createDimension('bins_across_track', nact)     # across track axis
    output_l2.createDimension('number_layers', nlay)         # layer axis

    # dimensions of two types
    _dims = ('bins_along_track', 'bins_across_track')
    _dims3d = ('bins_along_track', 'bins_across_track', 'number_layers',)

    # layer height
    _ = writevariablefromname(output_l2, 'layerheight', ('number_layers'),  data.zlay)
    # latitude
    _ = writevariablefromname(output_l2, 'latitude', _dims, data.lat)
    # longitude
    _ = writevariablefromname(output_l2, 'longitude', _dims, data.lon)

    # depending on gas write data
    if gas == "co2":
        xgas = "XCO2"
    elif gas == "ch4":
        xgas = "XCH4"
    elif gas == "no2":
        xgas = "XNO2"
    writegasandproxydata(xgas, output_l2, _dims, _dims3d, data.lvl2data, data.lvl2precision, data.avgkernel)
    output_l2.close()


def generate(sgmatmosphere_inputfile, level2_outputfile, gas, precision, seed=10):
    """Generate synthetic Level 2 data.

    Read SGM file (actual data) and create a synthetic level 2 data
    based on precision.

    Parameters
    ----------
    sgmatmosphere_inputfile : String
        SGM input file name.
    level2_outputfile : String
        Output file name of the Level 2 data.
    gas : String
        Kind of gas. Should be in lower-case.
    precision : Float
        Sigma value or level 2 precision. 0.2% should be given as
        0.002.
    seed : Int, Optional
        A seed to generate random numbers for level 2 precision.

    """
    # read data for a given gas
    data = readsgmatmosphere(sgmatmosphere_inputfile, gas)

    # generate level 2 data
    data.__setattr__("lvl2precision", data.actual_column*precision)
    data.__setattr__("lvl2data", getlevel2retrieval(data.actual_column, data.lvl2precision, seed))
    data.__setattr__("avgkernel", ones_like(data.columndensity))  # is set to one

    # save the level 2 data
    writedata_l2data(level2_outputfile, data, gas)
