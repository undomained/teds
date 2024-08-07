"""
Read-out Noise
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    shape = (cfg['dimensions'][dim_spat], cfg['dimensions'][dim_spec])
    gen.dtype = 'f8'
    gen.attr_names =  ['comment', 'long_name', 'units']
    gen.attr_vals = ['Calculation: sigma = sqrt(noise_g * signal + noise_n),\n Placeholder: random noise generated from single value parameter',
                     'read-out noise', 
                     'counts^2/e']

    # Calculate sigma
    ro_noise = cfg['read_out_noise']
    N0 = np.random.uniform(-ro_noise, ro_noise, size = shape)
    N0_scaled = N0**2
    gen.data = N0_scaled
        
    return gen
