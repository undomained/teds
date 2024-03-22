"""
Across track stage rotation
"""


import numpy as np
from tango.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    gen.attr_names =  ['long_name', 'units']
    gen.attr_vals = ['across track stage rotation per spatial sample', 'radians']
    
    dim_spat = 'spatial_samples_per_image'
    gen.dim_names = [dim_spat]
    shape = (cfg['dimensions'][dim_spat])
    gen.dtype = 'f8'
    
    fov = cfg['fov']
    fov_rad = fov*np.pi/180
    gen.data = np.linspace(-fov_rad/2, fov_rad/2, shape)

    return gen
