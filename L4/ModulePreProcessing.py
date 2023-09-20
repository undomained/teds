#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:54:22 2023

@author: manugv
"""

from .ModuleReadMicroHH import read_simulated_variable, get_interpolate_uv
from .ModuleLevel2 import readlevel2retrieval, generate_level2data
from netCDF4 import Dataset
from .ModuleDataContainer import DataCont
from pathlib import Path


mm_air = 0.0289647  # kg/mole
tot_air_moles = (101325.0 / 9.81) / mm_air  # moles
mm_co2 = 0.04401    # kg/mole
mm_ch4 = 0.01604    # kg/mole
mm_no2 = 0.0460055  # kg/mole


def readsgmatmosphere(sgmatmosphere_file, gas):
    f = Dataset(sgmatmosphere_file, "r")
    if gas == "co2":
        gascol = f["col_co2"][:].data
        dgascol = f["dcol_co2"][:].data
    elif gas == "ch4":
        gascol = f["col_ch4"][:].data
        dgascol = f["dcol_ch4"][:].data
    air = f["col_air"][:].data
    grid = DataCont()
    for ky in ["latitude", "longitude"]:
        grid.__setattr__(ky, f[ky][:])
    f.close()
    data = DataCont()
    data.__setattr__("actual_column", gascol)
    data.__setattr__("actual_column_air", air)
    data.__setattr__("dactual_column", dgascol)
    data.__setattr__("grid", grid)
    return data


def check_inputsimulationdata(param):
    # check for LSQ
    _file = param["simulationdata"]["sgm_atmosphere"]
    fl = Path(_file)
    if param["method"] == "LSQ":
        if not fl.is_file():
            print("Actual data is required for LSQ, please input a valid sgm_atmosphere")
            exit()

    # check level 2 data
    level2data = param["simulationdata"]["level2data"]
    if (level2data["generate"]):
        # if the data needs to be generated
        if not ((type(level2data["precision"]) == int) or (type(level2data["precision"]) == float)):
            print("Level 2 data cannot be generated as the 'precision' is not a number")
            exit()
        if not fl.is_file():
            print("Level 2 data cannot be generated as actual data is not present, please input a valid sgm_atmosphere")
            exit()
    else:
        # check if the files exist
        lvl2 = Path(level2data["file"])
        if not lvl2.is_file():
            print("Level 2 input file doesn't exist")
            exit()


def readyamlfile(filename):
    import yaml
    with open(filename, 'r') as file:
        param = yaml.safe_load(file)
    # check if simulation data files are proper
    check_inputsimulationdata(param)
    return param


def getendtoendsimdata(params):
    simparam = params["simulationdata"]
    level2data = params["simulationdata"]["level2data"]
    data = DataCont()
    # Read Actual data if exists
    if ~bool(simparam["sgm_atmosphere"]):
        data = readsgmatmosphere(simparam["sgm_atmosphere"], simparam["gas"])

    # level 2 data
    if (level2data["generate"]):
        generate_level2data(data, params["simulationdata"]["level2data"]["precision"])
    else:
        data = readlevel2retrieval(simparam["level2data"], simparam["gas"])
        data.grid.__setattr__("source", simparam["lat_lon_src"])

    # compute ppm to kg
    if simparam["gas"] == "co2":
        data.__setattr__("ppm_to_kg_gas", mm_co2*tot_air_moles/1e6)
    elif simparam["gas"] == "ch4":
        data.__setattr__("ppm_to_kg_gas", mm_ch4*tot_air_moles/1e6)
    elif simparam["gas"] == "no2":
        data.__setattr__("ppm_to_kg_gas", mm_no2*tot_air_moles/1e6)
    return data


def microhhvelocityinterp(params):
    prm = params["CFM"]["microhh"]
    data = read_simulated_variable(prm["path"], ['co2_m', 'u', 'v'], prm["time_stamp"])
    ht = params["CFM"]["plumeheight"]
    interpu, interpv = get_interpolate_uv(data, ht)
    return interpu, interpv, data
