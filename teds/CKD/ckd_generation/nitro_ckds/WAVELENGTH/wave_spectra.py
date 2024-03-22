"""
???
"""

import numpy as np
from tango.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()

    dim_spat = 'spatial_samples_per_image'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    N_spat_px = cfg['dimensions'][dim_spat]
    N_spec_px = cfg['dimensions'][dim_spec]
    gen.dtype = 'float32'

    gen.data = np.zeros((N_spat_px, N_spec_px))
    gen.attr_names.append('comment')
    gen.attr_vals.append('unused')

    return gen


