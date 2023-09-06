#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:54:22 2022

@author: manugv
"""
import os
import numpy as np
from netCDF4 import Dataset, num2date
import configparser
from scipy.interpolate import RegularGridInterpolator
from .ModuleDataContainer import DataCont


def _get_default_files(dir_data):
    """Gets all the .nc files in the folder

    Parameters
    ----------
    dir_data : string
        Directory of the data

    Returns
    -------
    List of strings
        List of file names
    """
    default_files = []
    for file in os.listdir(dir_data):
        if "default" in file:
            default_files.append(file)
    default_files.sort()
    return default_files


def _read_density(dir_data, _files, time):
    """Read density at a given time

    Parameters
    ----------
    dir_data : String
        Directory of the data
    _files : List of Strings
        List of file names with extension .nc
    time : Int64
        Time of the simulation of MicroHH

    Returns
    -------
    rho : Array
        Density at the given time
    """
    for _fl in _files:
        ff = Dataset(dir_data + _fl, "r")
        # id of time
        ix = np.where(np.isclose(ff["time"][:].data, time))[0]
        if len(ix) > 0:
            # rho for that time
            rho = ff["thermo/rho"][ix[0], :].data
            ff.close()
            return rho
        else:
            ff.close()
            continue


def _read_datatime(dir_data, _files, time):
    """Read density at a given time

    Parameters
    ----------
    dir_data : String
        Directory of the data
    _files : List of Strings
        List of file names with extension .nc
    time : Int64
        Time of the simulation of MicroHH

    Returns
    -------
    rho : Array
        Density at the given time
    """
    for _fl in _files:
        ff = Dataset(dir_data + _fl, "r")
        # id of time
        ix = np.where(np.isclose(ff["time"][:].data, time))[0]
        if len(ix) > 0:
            datetime = num2date(ff["time"][ix[0]], ff["time"].units)
            # rho for that time
            ff.close()
            return datetime
        else:
            ff.close()
            continue


def _get_ini_file(dir_data):
    """get ini file

    Parameters
    ----------
    dir_data : String
        Directory of the data

    Returns
    -------
    String
        File name of ini file
    """
    # find the file with ini extension in the folder
    for file in os.listdir(dir_data):
        if file.endswith(".ini"):
            return dir_data + file


def _get_domain(flname_ini):
    """Gets the grid of the given domain

    Parameters
    ----------
    flname_ini : string
        File name of the ini file

    Returns
    -------
    grid : data class
        Data class containing the grid of the microHH

    """

    # Create a grid container
    grid = DataCont()
    config = configparser.ConfigParser()
    config.read(flname_ini)
    grid.__setattr__("nz", config.getint("grid", "ktot"))
    grid.__setattr__("ny", config.getint("grid", "jtot"))
    grid.__setattr__("nx", config.getint("grid", "itot"))
    grid.__setattr__("zsize", config.getfloat("grid", "zsize"))
    grid.__setattr__("ysize", config.getfloat("grid", "ysize"))
    grid.__setattr__("xsize", config.getfloat("grid", "xsize"))
    grid.__setattr__("dz", grid.zsize / grid.nz)
    grid.__setattr__("dy", grid.ysize / grid.ny)
    grid.__setattr__("dx", grid.xsize / grid.nx)
    grid.__setattr__("x_nodes", np.linspace(0, grid.xsize, grid.nx + 1))
    grid.__setattr__("y_nodes", np.linspace(0, grid.ysize, grid.ny + 1))
    grid.__setattr__("z_nodes", np.linspace(0, grid.zsize, grid.nz + 1))
    grid.__setattr__("xc", 0.5 * (grid.x_nodes[1:] + grid.x_nodes[:-1]))
    grid.__setattr__("yc", 0.5 * (grid.y_nodes[1:] + grid.y_nodes[:-1]))
    grid.__setattr__("zc", 0.5 * (grid.z_nodes[1:] + grid.z_nodes[:-1]))
    return grid


def _get_source_strength(flname_ini, prefix):
    """get source and strength of the source

    Parameters
    ----------
    flname_ini : string
        Ini file
    prefix : string
        Type of gas

    Returns
    -------
    source : Vector[3]
        Source location [x,y,z]
    strength : Float64
        Strength of the emission in kilomoles/s
    """
    config = configparser.ConfigParser()
    config.read(flname_ini)
    # get source section and the index of the source
    source_sec = config["source"]
    srcs = source_sec.get("sourcelist")
    idx = (srcs.split(",")).index(prefix)

    # get source locations
    xloc = float((source_sec["source_x0"]).split(",")[idx])
    yloc = float((source_sec["source_y0"]).split(",")[idx])
    zloc = float((source_sec["source_z0"]).split(",")[idx])
    source = [xloc, yloc, zloc]
    # get strength
    strength = float((source_sec["strength"]).split(",")[idx])
    return source, strength


def _read_var(dir_data, prefix, time):
    """Read the data based on prefix

    Parameters
    ----------
    dir_data : String
        directory of the data
    prefix : String
        Gas that needs to be read.
    time : Int
        Time of the simulation

    Returns
    -------
    bdata : 1d-Array of Float64
        Array containing the binary data
    """
    with open(dir_data + prefix + "." + str(time).zfill(7)) as f:
        rectype = np.dtype(np.float64)
        bdata = np.fromfile(f, dtype=rectype)
        f.close()
    return bdata


def _get_variable(prefix, dir_data, ini_filename, time, dim):
    """Get source, strength and data from prefix at a given time.

    Parameters
    ----------
    prefix : String
        Gas that needs to be read.
    dir_data : String
        directory of the data
    ini_filename : String
        Ini file
    time : Int64
        Time of simulation in seconds
    dim : Array [3] Int64
        [nz, ny, nx] dimensions to reshape the binary data

    Returns
    -------
    Data : DataCont Class
        Data container with data
    """
    tmp = DataCont()
    # Get source from the ini file
    source, strength = _get_source_strength(ini_filename, prefix)
    tmp.__setattr__("source", source)
    tmp.__setattr__("strength", strength)
    bdata = _read_var(dir_data, prefix, time)
    tmp.__setattr__("conc", bdata.reshape(dim))
    return tmp


def read_simulated_variable(dir_data, prefix, time):
    """read a gas variable of microHH

    Parameters
    ----------
    dir_data : String
        directory of the data
    prefix : List of Strings or a String
        Gases/gas that need to be read.
    time : Int64
        Time of the simulation
        Example: co2.0034000 has prefix as "co2" and time as 340000
    Returns
    -------
    Data : class
        Container containing different variables.
    """
    data = DataCont()
    # Get the domain from the ini file
    ini_filename = _get_ini_file(dir_data)
    data.__setattr__("grid", _get_domain(ini_filename))

    # Read data
    dim = np.array([data.grid.nz, data.grid.ny, data.grid.nx])  # dimensions of binary data
    # Check if prefix is a list or one variable
    if isinstance(prefix, list):
        # Get data for each prefix
        for each_prefix in prefix:
            if (each_prefix == "u") or (each_prefix == "v"):
                veldata = _read_var(dir_data, each_prefix, time)
                data.__setattr__(each_prefix, veldata.reshape(dim))
            else:
                data.__setattr__(each_prefix, _get_variable(each_prefix, dir_data, ini_filename, time, dim))
    else:
        if (prefix == "u") or (prefix == "v"):
            veldata = _read_var(dir_data, prefix, time)
            data.__setattr__(prefix, veldata.reshape(dim))
        else:
            data.__setattr__(prefix, _get_variable(prefix, dir_data, ini_filename, time, dim))

    # get default files
    default_files = _get_default_files(dir_data)
    # Get time of the simulation
    data.__setattr__("datetime", _read_datatime(dir_data, default_files, time))
    # read density : rho
    data.__setattr__("density", _read_density(dir_data, default_files, time))
    return data


def get_interpolate_uv(data, plumeheight):
    """Get interpolation functions based on plumeheight in microHHdata


    Parameters
    ----------
    data : MicroHHdata class
        Class containing microHH data
    plumeheight : Float
        Height of the plume

    """
    iz = np.searchsorted(data.grid.zc, plumeheight)
    idz = data.grid.zc[iz] - plumeheight
    u1 = data.u[iz-1, :, :]*idz + data.u[iz, :, :]*(1-idz)
    v1 = data.v[iz-1, :, :]*idz + data.v[iz, :, :]*(1-idz)
    x1 = data.grid.xc - data.co2_m.source[0]
    y1 = data.grid.yc - data.co2_m.source[1]
    interp_u = RegularGridInterpolator((x1, y1), u1.T, bounds_error=False, fill_value=0)
    interp_v = RegularGridInterpolator((x1, y1), v1.T, bounds_error=False, fill_value=0)
    return interp_u, interp_v
