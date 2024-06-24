"""
Read-out Noise
"""

import numpy as np


def generate(ncc):
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_col']
    attr = {
        'long_name': 'read noise',
        'description' :  'noise model: sigma = sqrt(g*signal + n^2)', 
        'units': 'counts^2 e-1',
    }
    
    # Calculate sigma
    ro_noise = cfg['read_out_noise']
    N0 = np.random.uniform(-ro_noise, ro_noise, size = ncc.get_shape(dims))
    N0_scaled = N0**2
    ncc.create_var_auto(dims, N0_scaled, attr, 'f8')
