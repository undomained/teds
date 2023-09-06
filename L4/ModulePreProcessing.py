#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:54:22 2023

@author: manugv
"""

from .ModuleReadMicroHH import read_simulated_variable, get_interpolate_uv
from netCDF4 import Dataset
from ..lib.libNumTools import TransformCoords
from .ModuleDataContainer import DataCont


mm_air = 0.0289647  # kg/mole
tot_air_moles = (101325.0 / 9.81) / mm_air  # moles
mm_co2 = .04401    # kg/mole


def readsgmatmosphere(sgmatmosphere_file, gas):
    f = Dataset(sgmatmosphere_file, "r")
    if gas == "co2":
        gascol = f["col_co2"][:].data
        dgascol = f["dcol_co2"][:].data
    elif gas == "ch4":
        gascol = f["col_ch4"][:].data
        dgascol = f["dcol_ch4"][:].data
    air = f["col_air"][:].data
    f.close()
    return gascol, dgascol, air


def readlevel2retrieval(lvl2_file, gas):
    f = Dataset(lvl2_file, "r")
    data = DataCont()
    grid = DataCont()
    for ky in ["latitude", "longitude"]:
        grid.__setattr__(ky, f[ky][:])
    if gas == "co2":
        data.__setattr__("lvl2data", f["XCO2"][:])
        data.__setattr__("lvl2precision", f["precision XCO2"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XCO2"][:])
    elif gas == "no2":
        data.__setattr__("lvl2data", f["XNO2"][:])
        data.__setattr__("lvl2precision", f["precision XNO2"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XNO2"][:])
    elif gas == "ch4":
        data.__setattr__("lvl2data", f["XCH4"][:])
        data.__setattr__("lvl2precision", f["precision XCH4"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XCH4"][:])
    f.close()
    data.__setattr__("grid", grid)
    return data


def readyamlfile(filename):
    import yaml
    with open(filename, 'r') as file:
        param = yaml.safe_load(file)
    return param


def getendtoendsimdata(params):
    simparam = params["simulationdata"]
    data = readlevel2retrieval(simparam["level2data"], simparam["gas"])
    transform = TransformCoords(simparam["lat_lon_src"])
    xm, ym = transform.latlon2xymts(data.grid.latitude, data.grid.longitude)
    data.grid.__setattr__("xc", xm)
    data.grid.__setattr__("yc", ym)
    data.grid.__setattr__("source", simparam["lat_lon_src"])
    data.grid.__setattr__("sourcexy", [0, 0])
    # convert data to nodes from pixel centers
    # data.grid.__setattr__("x_nodes", xm)
    # data.grid.__setattr__("y_nodes", ym)
    data.__setattr__("source", simparam["lat_lon_src"])
    if ~bool(simparam["sgm_atmosphere"]):
        gas, dgas, air = readsgmatmosphere(simparam["sgm_atmosphere"], simparam["gas"])
        data.__setattr__("actual_column", gas)
        data.__setattr__("actual_column_air", air)
        data.__setattr__("dactual_column", dgas)

    data.__setattr__("ppm_to_kg_gas", mm_co2*tot_air_moles/1e6)
    return data


def microhhvelocityinterp(params):
    prm = params["CFM"]["microhh"]
    data = read_simulated_variable(prm["path"],
                                   ['co2_m', 'u', 'v'],
                                   prm["time_stamp"])
    ht = params["CFM"]["plumeheight"]
    interpu, interpv = get_interpolate_uv(data, ht)
    return interpu, interpv, data
