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