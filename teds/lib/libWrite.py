#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 2023.

@author: Manu Goudar
"""

from importlib.resources import files
from yaml import safe_load


consts_file = files("teds.lib").joinpath("constants_outputvariables.yaml")
with open(consts_file, "r") as file:
    variable_dict = safe_load(file)

# def writevariable(grp, data, _dims, _name, _long_name, _units="", _validmin=0, _validmax=0, _fillvalue=-32767):
#     gm_sza = grp.createVariable(_name, data.dtype, _dims)
#     gm_sza.long_name = _long_name
#     gm_sza.units = _units
#     gm_sza.valid_min = _validmin
#     gm_sza.valid_max = _validmax
#     gm_sza.FillValue = _fillvalue
#     gm_sza[:] = data


def writevariablefromname(grp, _name, dims, data):
    attr = variable_dict.get(_name)
    var = grp.createVariable(attr["name"], data.dtype, dims)
    var.long_name = attr["long_name"]
    var.units = attr["units"]
    var.valid_min = attr["valid_min"]
    var.valid_max = attr["valid_max"]
    var.FillValue = attr["FillValue"]
    var[:] = data
    return var
