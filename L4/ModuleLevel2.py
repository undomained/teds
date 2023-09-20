#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:54:22 2023

@author: manugv
"""

from netCDF4 import Dataset
from .ModuleDataContainer import DataCont
from numpy import random


def readlevel2retrieval(lvl2_file, gas, data):
    f = Dataset(lvl2_file, "r")
    # if grid attribute doesn't exist
    if not hasattr(data, "grid"):
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
    # add grid
    data.__setattr__("grid", grid)
    return data


def getlevel2retrieval(conc, precision, seed=10):
    noise = random.default_rng(seed).normal(0, 1, conc.shape)
    return conc + noise*precision


def generate_level2data(data, precision):
    data.__setattr__("lvl2precision", data.actual_column*precision)
    data.__setattr__("lvl2data", getlevel2retrieval(data.actual_column, data.lvl2precision))
    data.__setattr__("avg_kernel", 1.0)
