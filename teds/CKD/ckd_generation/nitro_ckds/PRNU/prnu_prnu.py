"""
Pixel Response Non-Uniformity
Does not include Quantum Efficiency
"""

import numpy as np
from tango.ckd_generation.generator_class import *
from scipy.ndimage import gaussian_filter

def generate(cfg):
    gen = ckd_generator()
    dim_spat = 'spatial_detector_pixels'
    dim_spec = 'spectral_detector_pixels'
    gen.dim_names = [dim_spat, dim_spec]
    N_spat_px = cfg['dimensions'][dim_spat]
    N_spec_px = cfg['dimensions'][dim_spec]
    gen.dtype = 'f8'

    prnu_budget = cfg['prnu'] # in %, from TANGO Nitro SNR budget (TNO)

    # Create a T image, temporary of course, until we get new data
    patternsize = 120
    barcols = int(patternsize/2)
    barnegcols = patternsize - barcols
    barrow = np.roll(np.append(np.ones((barcols)), np.zeros((barnegcols))), int(barnegcols/2))
    
    weight = int(patternsize/10)
    pillarnegcols = patternsize - weight
    pillarrow  = np.roll(np.append(np.ones((weight)), np.zeros((pillarnegcols))), int(pillarnegcols/2))

    pattern = np.zeros((patternsize, patternsize))
    marginrows = 20
    for r in range(patternsize)[marginrows:-marginrows]:
        if r >= patternsize - (marginrows + weight):
            pattern[r] = barrow
        else:
            pattern[r] = pillarrow
    pattern = 1 - pattern # inverse


    result = np.zeros((N_spat_px, N_spec_px), dtype=pattern.dtype)
    for i in range(N_spat_px):
        for j in range(N_spec_px):
            result[i, j] = pattern[i % patternsize, j % patternsize]
    
    result = gaussian_filter(result, 8, order=0, output=None)
    noise = np.random.rand(N_spat_px, N_spec_px)
    result += noise

    # Scale
    prnu_max = 1 + prnu_budget/100
    prnu_min = 1 - prnu_budget/100
    result = result-np.min(result)
    result = result /  (np.max(result)) * (prnu_max - prnu_min)
    result = result + prnu_min
    gen.data = result

    gen.attr_names =  ['comment', 'long_name', 'units']
    gen.attr_vals = ['Placeholder, does not include Quantum Efficiency',
                     'pixel response non-uniformity', 
                     '1']

    return gen
