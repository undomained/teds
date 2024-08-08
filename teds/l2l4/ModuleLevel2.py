#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:54:22 2023

@author: manugv
"""

from netCDF4 import Dataset
from .ModuleDataContainer import DataCont
import sys

def read_level2_product(lvl2_file, gas):

    data = DataCont()
    f = Dataset(lvl2_file, "r")
    
    data.__setattr__("lat", f["lat"][:])
    data.__setattr__("lon", f["lon"][:])
    if gas == "co2":
        data.__setattr__("Xgas", f["XCO2 proxy"][:])
        data.__setattr__("Xgas_precision", f["precision XCO2 proxy"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XCO2"][:])
    elif gas == "no2":
        data.__setattr__("Xgas", f["XNO2"][:])
        data.__setattr__("Xgas_precision", f["precision XNO2"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XNO2"][:])
    elif gas == "ch4":
        data.__setattr__("Xgas", f["XCH4"][:])
        data.__setattr__("Xgas_precision", f["precision CH4"][:])
        data.__setattr__("avg_kernel", f["col avg kernel XCH4"][:])
    f.close()
  
    return data
