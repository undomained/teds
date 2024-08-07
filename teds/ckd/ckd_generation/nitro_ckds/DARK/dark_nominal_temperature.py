"""
Dark Current
"""

import numpy as np
from teds.CKD.ckd_generation.generator_class import *

def generate(cfg):
    gen = ckd_generator()
    gen.attr_names =  ['long_name', 'units', 'comment']
    gen.attr_vals = ['temperature', 'degr (C)', 'Placeholder, no data available']
    
    dim= 'single_double'
    gen.dim_names = [dim]
    gen.dtype = 'f8'
    dark_nominal_temperature = 15
    gen.data = dark_nominal_temperature

    return gen
