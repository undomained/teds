"""
Read-out Noise
"""

import numpy as np


def generate(ncc):
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']
    attr = {
        'long_name': 'n2: read noise',
        'description' :  'noise model: sigma = sqrt(g*signal + n^2)', 
        'units': 'counts^2 e-1',
    }
    
    # Calculate sigma
    ro_noise = cfg['read_out_noise']
    N0 = np.ones(ncc.get_shape(dims)) * ro_noise
    N0_scaled = N0**2
    ncc.create_var_auto(dims, N0_scaled, attr, 'f8')
