"""
Dark Current
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    gen.attr_names =  ['long_name', 'units', 'comment']
    gen.attr_vals = ['detector offset', 'counts', 'Placeholder, no data available']
    
    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    shape = (cfg['dimensions'][dim_spat], cfg['dimensions'][dim_spec])
    gen.dtype = 'f8'
    
    dark_offset = np.zeros(shape)
    gen.data = dark_offset

    return gen
