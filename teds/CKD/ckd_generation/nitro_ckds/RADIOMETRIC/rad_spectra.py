"""
Optical Radiometric Response
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    gen.attr_names =  ['long_name', 'units', 'comment']
    gen.attr_vals = ['optical radiometric response', '1', 
                     'Placeholder, single value matrix']
    
    gen.dim_names = ['spatial_samples_per_image',  'spectral_detector_pixels']
    shape = [cfg['dimensions'][dimname] for dimname in gen.dim_names]
    gen.dtype = 'f8'
    
    c_rad = cfg['c_opt_rad_resp']
    C = np.ones(shape) * c_rad
    gen.data = C

    return gen
