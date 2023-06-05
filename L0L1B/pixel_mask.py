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

def determine_pixel_mask(paths, filen_pixel_mask, file1b, file1b_ref):
    
    
    path = paths.project+paths.data_interface+paths.interface_l1b

    input = nc.Dataset(path+file1b, mode='r')
    radiance= deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    input.close()
    input = nc.Dataset(path+file1b_ref, mode='r')
    radiance_ref= deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    input.close()

    diff = (radiance[0,:,:]-radiance_ref[0,:,:])/radiance_ref[0,:,:]*100.
    mask = (np.abs(diff) < 3.).data
    
    # l1b_data['noise'] = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])
    np.save(path+filen_pixel_mask, mask)

    return