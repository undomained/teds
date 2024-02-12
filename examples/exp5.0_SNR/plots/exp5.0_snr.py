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

nsza  = 1
nalb  = 41
nact  = 100
nwave = 550

szas = [70.]
albs = np.arange(nalb)*0.01 + 0.01
snr_max= np.zeros((nalb,nsza))
snr_min= np.zeros((nalb,nsza))
rad_max= np.zeros((nalb,nsza))
rad_min= np.zeros((nalb,nsza))
l2prec_median = np.zeros((nalb,nsza))
co2 = np.zeros((nact, nsza))

dum1 = np.zeros(nact)
dum2 = np.zeros(nact)
for ialb, alb in enumerate(albs):
    for isza, sza in enumerate(szas):
    
        filenamel1b = path_interface + 'level1b/Tango_Carbon_l1b_exp5.0_sza'+"{:.1f}".format(sza)+'_alb'+"{:.3f}".format(alb)+'.nc'
        filenamel2  = path_interface + 'level2/Tango_Carbon_l2_exp5.0_sza'+"{:.1f}".format(sza) +'_alb'+"{:.3f}".format(alb)+'.nc'

        l1b_data = nc.Dataset(filenamel1b)
        l2_data  = nc.Dataset(filenamel2)

        wave  = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['wavelength'][:]))
        rad   = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['radiance'][:]))
        noise = np.array(deepcopy(l1b_data['OBSERVATION_DATA']['radiance_noise'][:]))
        snr = rad/noise
    
        lmax = 114
        lmin = 502
        snr_max[ialb, isza] = np.median(snr[0,:,lmax]) 
        snr_min[ialb, isza] = np.median(snr[0,:,lmin])
        rad_max[ialb, isza] = np.median(rad[0,:,lmax]) 
        rad_min[ialb, isza] = np.median(rad[0,:,lmin])
        
        prec_xco2_proxy = np.array(deepcopy(l2_data['precision XCO2 proxy'][:]).flatten())

        l2prec_median[ialb,isza] = np.median(prec_xco2_proxy)

        # if np.abs(alb-0.15) < 1.e-5:
            
        #     fig = plt.figure(figsize=(26, 7), dpi=100,)
        #     ax1 = fig.add_subplot(121)
            
        #     ax1.plot(wave[0,:],rad[0,50,:])
        #     ax1.set_xlim([1590,1685])
        #     ax1.set_xlabel('wavelength [nm]')
        #     ax1.set_ylabel('radiance [ph./(s sr m$^2$ nm)]')

        #     ax2 = fig.add_subplot(122)
            
        #     ax2.plot(wave[0,:],snr[0,50,:])
        #     ax2.set_xlim([1590,1685])
        #     ax2.set_xlabel('wavelength [nm]')
        #     ax2.set_ylabel('snr [1]')

        #     plt.savefig('plots/radiance_snr.png',)
        #     sys.exit()

l2_req   = 3.8
wave_max = wave[0,lmax] 
wave_min = wave[0,lmin]
req_snr_max = np.interp(l2_req, np.flip(l2prec_median[:,0]), np.flip(snr_max[:, isza]))
req_snr_min = np.interp(l2_req, np.flip(l2prec_median[:,0]), np.flip(snr_min[:, isza]))
req_rad_max = rad_max[14, isza]
req_rad_min = rad_min[14, isza]

sw_fig1 = False
sw_fig2 = True

if(sw_fig1):
    fig = plt.figure(figsize=(26, 7), dpi=100,)
    ax1 = fig.add_subplot(121)
    ax1.plot(albs/0.15,snr_max[:, isza], color = 'green', label = 'SNR$_\mathrm{max}$')
    ax1.plot(albs/0.15,snr_min[:, isza], color = 'blue', label = 'SNR$_\mathrm{min}$')
    ax1.set_xlabel('scaling factor $s$ [1]')
    ax1.set_ylabel('SNR [1]')
    plt.legend()
    
    sys.exit()
    req_scaling = np.interp(l2_req, np.flip(l2prec_median[:,0]), np.flip(albs))/0.15
    print('required scaling factor:', req_scaling)
    
    ax2 = fig.add_subplot(122)
    ax2.plot(albs/0.15,l2prec_median[:, isza], color = 'grey')
    ax2.set_xlabel('scaling factor $s$ [1]')
    ax2.set_ylabel('XCO$_2^\mathrm{prox}$ [ppm]')
    
    textstr = '\n'.join((
        r'required L2 precsion =%.2f' % (l2_req, ) + ' ppm',
        r'degradation factor $s =%.2f$' % (req_scaling, )))
    props = dict(boxstyle='round', facecolor='grey', alpha=0.2)
    ax2.text(0.15, 0.95, textstr, transform=ax2.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    plt.savefig('plots/l2_precision.png',)
    
if(sw_fig2):

    fig = plt.figure(figsize=(11, 7), dpi=100,)
    
    textstr_max = '\n'.join((
        r'$\lambda_\mathrm{max}=%.2f$' % (wave_max, ) + 'nm',
        r'$I_\mathrm{max}=%.2E$' % (req_rad_max, ) + ' ph/(s sr m$^2$ nm)',
        r'$SNR_\mathrm{max}=%.2f$' % (req_snr_max, )))
    
    textstr_min = '\n'.join((
        r'$\lambda_\mathrm{min}=%.2f$' % (wave_min, ) + 'nm',
        r'$I_\mathrm{min}=%.2E$' % (req_rad_min, ) + ' ph/(s sr m$^2$ nm)',
        r'$SNR_\mathrm{min}=%.2f$' % (req_snr_min, )))
    
    ax1 = fig.add_subplot(111)
    ax1.plot(l2prec_median[:,0],snr_max[:, isza], color = 'green', label = 'SNR$_\mathrm{max}$')
    ax1.plot(l2prec_median[:,0],snr_min[:, isza], color = 'blue', label = 'SNR$_\mathrm{min}$')
    ax1.set_xlabel('XCO$_2^\mathrm{prox}$ [ppm]')
    ax1.set_ylabel('SNR [1]')
    plt.legend()
    props = dict(boxstyle='round', facecolor='green', alpha=0.2)
    ax1.text(0.15, 0.95, textstr_max, transform=ax1.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    props = dict(boxstyle='round', facecolor='blue', alpha=0.2)
    ax1.text(0.15, 0.65, textstr_min, transform=ax1.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    plt.savefig('plots/snr.png',)
    
plt.show()
