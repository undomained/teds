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
path_interface ='/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/'

dimt1 = 21

co2_T10  = np.zeros((dimt1))
temperature_range1 = (np.arange(dimt1)*0.01 -0.1)
scale = 10.

for it, delta_temperature in enumerate(temperature_range1):
    
    target_temperature1 = delta_temperature + 10.
    print('===============================================================')
    print('temperature: ', f"{target_temperature1: 04.2f}")
    print('===============================================================')

    filename = path_interface + 'level2/Tango_Carbon_l2_exp4.1_'+\
                'temp_'+f"{target_temperature1:05.2f}" + \
                '_scale_' + f"{scale:04.1f}" + ".nc"

    l2_data = nc.Dataset(filename)
    xco2_proxy = deepcopy(l2_data['XCO2 proxy'][:]).flatten()
    l2_data.close()
    mask = np.logical_not(np.isnan(xco2_proxy))

    co2_T10[it] = np.mean(xco2_proxy[mask])    

dimt2 = 21
co2_T15  = np.zeros((dimt2))
temperature_range2 = (np.arange(dimt2)*0.01 -0.1)

for it, delta_temperature in enumerate(temperature_range2):
    target_temperature2 = delta_temperature + 15.
    print('===============================================================')
    print('temperature: ', f"{target_temperature2: 04.2f}")
    print('===============================================================')

    filename = path_interface + 'level2/Tango_Carbon_l2_exp4.1_'+\
                'temp_'+f"{target_temperature2:05.2f}" + \
                '_scale_' + f"{scale:04.1f}" + ".nc"

    l2_data = nc.Dataset(filename)
    xco2_proxy = deepcopy(l2_data['XCO2 proxy'][:]).flatten()
    l2_data.close()
    mask = np.logical_not(np.isnan(xco2_proxy))

    co2_T15[it] = np.mean(xco2_proxy[mask])
    
fig = plt.figure(figsize=(12, 7), dpi=100,)
ax1 = fig.add_subplot(111)
ax1.plot(temperature_range1*1.E3, co2_T10-co2_T10[10], marker='+', label='T = 10$^o$C', color = 'red')
ax1.plot(temperature_range2*1.E3, co2_T15-co2_T15[10], marker='+', label='T = 15$^o$C', color = 'blue')
ax1.set_title('Temperature induced XCO$_2$ error, Owl640-S low gain')
ax1.set_ylabel('XCO$_2$ proxy bias [ppm]')
ax1.set_xlabel('detector temperature error $\delta$T [mK]')
plt.legend()
plt.show()
plt.savefig('plots/exp4.1_det_temp.png',)


plt.show()
