"""
???
"""

import numpy as np
from tango.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()

    gen.dim_names = 'number_of_views'
    gen.dtype = 'u1'
    
    gen.attr_names = ['comment']
    gen.attr_vals = ['detector pixel mask']
    gen.attr_vals.append((np.uint8(0), np.uint8(1)))
    gen.attr_vals.append('good bad')
    gen.attr_vals.append('unused')    

    gen.data = []
    
    return gen