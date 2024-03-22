"""
Detector Pixel Mask
"""
import numpy as np
from tango.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()

    gen.dim_names = ['spatial_detector_pixels', 'spectral_detector_pixels']
    dim0 = cfg['dimensions'][gen.dim_names[0]]
    dim1 = cfg['dimensions'][gen.dim_names[1]]
    gen.dtype = 'u1'
    mask = np.zeros((dim0, dim1))
    
    gen.attr_names = ['long_name', 'flag_values', 'flag_meanings']
    gen.attr_vals = ['detector pixel mask']
    gen.attr_vals.append((np.uint8(0), np.uint8(1)))
    gen.attr_vals.append('good bad')

    gen.data = mask

    return gen