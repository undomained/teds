"""
Spatial Sample to Pixel Conversion
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    gen.attr_names =  ['long_name', 'units', 'comment']
    gen.attr_vals = ['floating point spatial detector pixel indices of the spectra', '1', 
                     'for conversion from spatial sample to spatial pixel and vice-versa']
    
    gen.dim_names = ['spatial_samples_per_image',  'spectral_detector_pixels']
    shape = [cfg['dimensions'][dimname] for dimname in gen.dim_names]
    gen.dtype = 'f8'
    
    Npx_spat = cfg['dimensions']['spatial_detector_pixels']
    row = np.linspace(Npx_spat, 0, shape[0]) # spatial pixel index for each spatial sample

    gen.data = np.repeat(np.array(row, ndmin=2), shape[1], axis=0).T

    return gen