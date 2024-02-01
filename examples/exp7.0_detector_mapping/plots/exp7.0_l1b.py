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
    
filenamel1b = path_interface + 'level1b/Tango_Carbon_l1b_exp7.0.nc'

l1b_data = nc.Dataset(filenamel1b)

wave  = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['wavelength'][:]))
rad   = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['radiance'][:]))
noise = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['radiance_noise'][:]))

l1b_data.close()

nact = 100
alb = np.zeros(nact) + 0.10
alb[::2] = 0.60

iact = 97
fig = plt.figure(figsize=(26, 7), dpi=100,)
ax1 = fig.add_subplot(121)

ax1.plot(wave[0,:],rad[0,iact-1,:]/rad[0,iact-1,114])
ax1.plot(wave[0,:],rad[0,iact,:]/rad[0,iact,114])
ax1.set_xlim([1590,1685])
ax1.set_xlabel('wavelength [nm]')
ax1.set_ylabel('radiance [ph./(s sr m$^2$ nm)]')

ax2 = fig.add_subplot(122)

ax2.plot(wave[0,:],rad[0,iact,:]/noise[0,iact,:])
ax2.set_xlim([1590,1685])
ax2.set_xlabel('wavelength [nm]')
ax2.set_ylabel('snr [1]')

plt.savefig('plots/spectra_iact'+"{:.0f}".format(iact) + '.png',)


