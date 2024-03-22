"""
Noise Conversion Gain
"""

import numpy as np
from tango.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    shape = (cfg['dimensions'][dim_spat], cfg['dimensions'][dim_spec])
    gen.dtype = 'f8'

    gain = cfg['noise_conversion_gain']
    gen.data = np.ones(shape) * gain

    gen.attr_names =  ['comment', 'long_name', 'units']
    gen.attr_vals = ['Placeholder, value based on TANGO Carbon ckd',
                     'noise conversion gain', 
                     '1']

    return gen
