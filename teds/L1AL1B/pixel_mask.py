#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 16:11:25 2023

@author: jochen
"""
import netCDF4 as nc
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
import sys
import random

def determine_pixel_mask(filen_pixel_mask, file1b, file1b_ref):
    
    print(file1b)
    print(file1b_ref)
    input = nc.Dataset(file1b, mode='r')
    radiance= deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    input.close()
    input = nc.Dataset(file1b_ref, mode='r')
    radiance_ref= deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    input.close()

    diff = (radiance[0,:,:]-radiance_ref[0,:,:])/radiance_ref[0,:,:]*100.
    mask = (np.abs(diff) < 3.).data
    
    # for iact in range(diff[:,0].size):
    #     print(iact, np.sum(mask[iact,:]))
    # fig = plt.figure(figsize=(12, 12), dpi=100)
    # m = plt.imshow(diff, vmin = -10, vmax = 10)
    # plt.title(np.sum(mask)/diff.size)
    # cbar = plt.colorbar(m, orientation = 'vertical', label = 'counts [bu]', fraction=0.03, pad=0.1)
    # sys.exit()
    np.save(filen_pixel_mask, mask)

    return


def determine_RTS_pixel_mask(filen_pixel_mask, file1b, nbin, seed, fraction):
    
    print(file1b)
    input = nc.Dataset(file1b, mode='r')
    radiance= deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    input.close()

    nalt = radiance[:,0,0].size
    nact = radiance[0,:,0].size
    nwav = radiance[0,0,:].size
    
    pixel_binned = np.zeros([nalt,nact,nwav])
    
    for ialt in range(nalt):
        pixel_unbinned = np.zeros([nbin*nact*nwav])
        
        random.seed(seed)  #define seed 
        print('fraction of random telegraph pixels: ', fraction*pixel_unbinned.size)
        indices = random.sample(range(pixel_unbinned.size), int(fraction * pixel_unbinned.size))
        # set value to 1 for the selected pixels
        pixel_unbinned[indices] = 1.
        pixel_unbinned = np.reshape(pixel_unbinned, (nact,nwav,nbin))
        pixel_binned[ialt,:,:] = np.sum(pixel_unbinned, axis =2)

    mask = (pixel_binned[0,:,:] < 1.).data
    
    print('Fraction of bad binned pixels: ', (1-np.sum(mask)/(nact*nwav))*100, ' %')
    
#        fig = plt.figure(figsize=(12, 12), dpi=100)
#        m = plt.imshow(mask, vmin = -0, vmax = 1, aspect = 5)
#        plt.title(np.sum(pixel_unbinned)/pixel_unbinned.size)
#        cbar = plt.colorbar(m, orientation = 'vertical', label = 'counts [bu]', fraction=0.03, pad=0.1)

    np.save(filen_pixel_mask, mask)
    return