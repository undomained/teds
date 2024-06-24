"""
Noise Conversion Gain
"""

import numpy as np


def generate(ncc):
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_col']
    gain = cfg['noise_conversion_gain']
    g = np.ones(ncc.get_shape(dims)) * gain
    
    attr = {
        'long_name': 'noise conversion gain',
        'description' :  'noise model: sigma = sqrt(g*signal + n^2)', 
        'units': 'counts e-1',
        'comment': 'Placeholder, value based on Tango Carbon ckd'
    }
    ncc.create_var_auto(dims, g, attr, 'f8')
