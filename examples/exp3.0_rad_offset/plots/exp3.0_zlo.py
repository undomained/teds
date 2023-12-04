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

nzlo = 11
nsza = 8
zlo = np.arange(nzlo)*0.005
sza = [70., 60, 50, 40, 30, 20, 10, 0] 
co2 = np.zeros((nzlo, nsza))
co2_ns = np.zeros((nzlo, nsza))
ch4_ns = np.zeros((nzlo, nsza))
for izlo, rad_offset in enumerate(zlo):
    filename = path_interface + 'level2/Tango_Carbon_l2_exp3.0_rad_offset'+ \
                "{:.3f}".format(rad_offset)+'.nc'    
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
    co2[izlo,:] = xco2_proxy[:]
    co2_ns[izlo,:] = xco2_ns[:]
    ch4_ns[izlo,:] = xch4_ns[:]

delta_xco2    = np.zeros((nzlo, nsza))
delta_xco2_ns = np.zeros((nzlo, nsza))
delta_xch4_ns = np.zeros((nzlo, nsza))

for izlo in range(nzlo):
    delta_xco2[izlo,:]=co2[izlo,:]-co2[0,:]
    delta_xco2_ns[izlo,:]=(co2_ns[izlo,:]-co2_ns[0,:])/co2_ns[0,:]*100.
    delta_xch4_ns[izlo,:]=(ch4_ns[izlo,:]-ch4_ns[0,:])/ch4_ns[0,:]*100.
    
fig = plt.figure(figsize=(26, 7), dpi=100,)
ax1 = fig.add_subplot(131)
mesh1 = ax1.pcolormesh(sza, zlo*100., delta_xco2_ns) #, vmin=390., vmax=435., cmap = 'BuPu', transform=ccrs.PlateCarree(), alpha = 0.8)
cbar1 = plt.colorbar(mesh1, ax=ax1, orientation='vertical', fraction=0.04, pad=0.05)
cbar1.set_label('$\delta$ XCO$_2$ [%]')
ax1.set_title('non-scattering XCO$_2$ error: ZLO')
ax1.set_ylabel('ZLO [%]')
ax1.set_xlabel('SZA [degree]')

ax1 = fig.add_subplot(132)
mesh1 = ax1.pcolormesh(sza, zlo*100., delta_xch4_ns) #, vmin=390., vmax=435., cmap = 'BuPu', transform=ccrs.PlateCarree(), alpha = 0.8)
cbar1 = plt.colorbar(mesh1, ax=ax1, orientation='vertical', fraction=0.04, pad=0.05)
cbar1.set_label('$\delta$ XCH$_4$ [%]')
ax1.set_title('non-scattering XCH$_4$ error: ZLO')
ax1.set_ylabel('ZLO [%]')
ax1.set_xlabel('SZA [degree]')

ax1 = fig.add_subplot(133)
mesh1 = ax1.pcolormesh(sza, zlo*100., delta_xco2) #, vmin=390., vmax=435., cmap = 'BuPu', transform=ccrs.PlateCarree(), alpha = 0.8)
cbar1 = plt.colorbar(mesh1, ax=ax1, orientation='vertical', fraction=0.04, pad=0.05)
cbar1.set_label('$\delta$ XCO$_2$ [ppm]')
ax1.set_title('proxy XCO$_2$ error: ZLO')
ax1.set_ylabel('ZLO [%]')
ax1.set_xlabel('SZA [degree]')

fig.tight_layout(pad = 5.)
plt.show()
plt.savefig('plots/zlo.png',)


plt.show()
