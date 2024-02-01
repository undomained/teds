#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 21:06:07 2023

@author: jochen
"""
import netCDF4 as nc
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 14,})
path_interface ='/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/'
    
filenamel2 = path_interface + 'level2/Tango_Carbon_l2_exp7.0.nc'
l2_data = nc.Dataset(filenamel2)
        
prec_xco2_proxy = np.array(deepcopy(l2_data['precision XCO2 proxy'][:]).flatten())
xco2_proxy = np.array(deepcopy(l2_data['XCO2 proxy'][:]).flatten())

l2_data.close()
plt.plot(xco2_proxy)
plt.show()
