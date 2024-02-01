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

dimt = 16
co2      = np.zeros((dimt))
prec_co2 = np.zeros((dimt))
co2_dc      = np.zeros((dimt))
prec_co2_dc = np.zeros((dimt))
co2_qe      = np.zeros((dimt))
prec_co2_qe = np.zeros((dimt))

temperature_range = np.arange(dimt)*1.0
scale = 10.

for it, target_temperature in enumerate(temperature_range):
     
    print('===============================================================')
    print('temperature: ', f"{target_temperature: 04.2f}")
    print('===============================================================')

    filename = path_interface + 'level2/Tango_Carbon_l2_exp4.0_'+\
                'temp_'+f"{target_temperature:05.2f}" + \
                '_scale_' + f"{scale:04.1f}" + ".nc"

    l2_data = nc.Dataset(filename)
    xco2_proxy = deepcopy(l2_data['XCO2 proxy'][:]).flatten()
    xch4_proxy = deepcopy(l2_data['XCH4 proxy'][:]).flatten()
    xco2_ns = deepcopy(l2_data['XCO2'][:]).flatten()
    xch4_ns = deepcopy(l2_data['XCH4'][:]).flatten()
    prec_xco2_proxy = deepcopy(l2_data['precision XCO2 proxy'][:]).flatten()
    prec_xch4_proxy = deepcopy(l2_data['precision XCH4 proxy'][:]).flatten()
    prec_xco2_ns = deepcopy(l2_data['precision XCO2'][:]).flatten()
    prec_xch4_ns = deepcopy(l2_data['precision XCH4'][:]).flatten()
    l2_data.close()
    mask = np.logical_not(np.isnan(xco2_proxy))

    co2[it] = np.mean(xco2_proxy[mask])
    prec_co2[it] = np.mean(prec_xco2_proxy[mask])
    
for it, target_temperature in enumerate(temperature_range):
     
    filename = path_interface + 'level2/Tango_Carbon_l2_exp4.0_'+\
                'temp_'+f"{target_temperature:05.2f}" + \
                '_scale_' + f"{scale:04.1f}" + "_dc.nc"

    l2_data = nc.Dataset(filename)
    xco2_proxy = deepcopy(l2_data['XCO2 proxy'][:]).flatten()
    xch4_proxy = deepcopy(l2_data['XCH4 proxy'][:]).flatten()
    xco2_ns = deepcopy(l2_data['XCO2'][:]).flatten()
    xch4_ns = deepcopy(l2_data['XCH4'][:]).flatten()
    prec_xco2_proxy = deepcopy(l2_data['precision XCO2 proxy'][:]).flatten()
    prec_xch4_proxy = deepcopy(l2_data['precision XCH4 proxy'][:]).flatten()
    prec_xco2_ns = deepcopy(l2_data['precision XCO2'][:]).flatten()
    prec_xch4_ns = deepcopy(l2_data['precision XCH4'][:]).flatten()
    l2_data.close()
    mask = np.logical_not(np.isnan(xco2_proxy))

    co2_dc[it] = np.mean(xco2_proxy[mask])
    prec_co2_dc[it] = np.mean(prec_xco2_proxy[mask])

for it, target_temperature in enumerate(temperature_range):
     
    filename = path_interface + 'level2/Tango_Carbon_l2_exp4.0_'+\
                'temp_'+f"{target_temperature:05.2f}" + \
                '_scale_' + f"{scale:04.1f}" + "_qe.nc"

    l2_data = nc.Dataset(filename)
    xco2_proxy = deepcopy(l2_data['XCO2 proxy'][:]).flatten()
    xch4_proxy = deepcopy(l2_data['XCH4 proxy'][:]).flatten()
    xco2_ns = deepcopy(l2_data['XCO2'][:]).flatten()
    xch4_ns = deepcopy(l2_data['XCH4'][:]).flatten()
    prec_xco2_proxy = deepcopy(l2_data['precision XCO2 proxy'][:]).flatten()
    prec_xch4_proxy = deepcopy(l2_data['precision XCH4 proxy'][:]).flatten()
    prec_xco2_ns = deepcopy(l2_data['precision XCO2'][:]).flatten()
    prec_xch4_ns = deepcopy(l2_data['precision XCH4'][:]).flatten()
    l2_data.close()
    mask = np.logical_not(np.isnan(xco2_proxy))

    co2_qe[it] = np.mean(xco2_proxy[mask])
    prec_co2_qe[it] = np.mean(prec_xco2_proxy[mask])

fig = plt.figure(figsize=(12, 7), dpi=100,)
ax1 = fig.add_subplot(111)
ax1.plot(temperature_range, prec_co2, marker='+', label='EoL', color = 'red')
ax1.plot(temperature_range, prec_co2_dc, marker='+', label='EoL, prnu*qe =1', color = 'blue')
ax1.plot(temperature_range, prec_co2_qe, marker='+', label='EoL, dc, doff =0', color = 'green')
ax1.set_title('XCO$_2$ sensitivity to detector temperature T, Owl640-S low gain')
ax1.set_ylabel('XCO$_2$ proxy precision [ppm]')
ax1.set_xlabel('detector temperature T [$^o$C]')
plt.legend()
plt.show()
plt.savefig('plots/exp4.0_det_temp_comp.png',)


plt.show()
