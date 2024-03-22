"""
Quantum Efficiency
Wavelength dependend -> Uses Wavelength Map
"""

import numpy as np
from scipy.interpolate import interp1d
from tango.ckd_generation.generator_class import *
import os

def generate(cfg):
    gen = ckd_generator()

    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    dim_names = [dim_spat, dim_spec]

    # data, source: TANGO Nitro SNR budget.xlsx (TNO)
    wl = [405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490]
    qe = [0.54, 0.56, 0.58, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.685, 0.70, 0.7, 0.705, 0.71, 0.71]
    f = interp1d(wl,qe, fill_value = 'extrapolate')

    # check if there is a generation script for wavelength map
    wl_map_path = cfg['paths']['dir_nitro'] + '/WAVELENGTH/wave_map.py'
    if os.path.isfile(wl_map_path):
        mod = import_mod(wl_map_path)
        wl_gen = mod.generate(cfg)
        wl_map = wl_gen.data
        QE = f(wl_map)
        gen.data = QE
    else:
        gen.data = np.ones((cfg['dimensions']['spatial_detector_pixels'], cfg['dimensions']['spectral_detector_pixels']))
    
    gen.attr_names = ['comment', 'units']
    gen.attr_vals  = ['Generated using wavelength map and QE data from TNO', '1']
    gen.dtype = 'float32'
    gen.dim_names = dim_names

    return gen






