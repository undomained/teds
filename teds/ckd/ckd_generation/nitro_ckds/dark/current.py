"""
Dark Current
"""

import numpy as np

def generate(ncc):

    cfg = ncc.cfg
    seed = cfg['seed']
    mu = cfg['dark_current']
    sigma is None
    if 'dark_current_sigma' in cfg:
        sigma = cfg['dark_current_sigma']

    attr = {
        'long_name': 'detector dark current',
        'units': 'counts/s',
        'comment': 'Not temperature dependent'
    }
    
    dims = ['detector_row', 'detector_column']
    shape = ncc.get_shape(dims)

    if sigma is not None:
        N = shape[0]*shape[1]
        if seed is not None:
            np.random.seed(seed)
    #    s = np.random.lognormal(mu, sigma, N)    
        s = np.random.normal(mu, sigma, N)    
        dark_current = s.reshape(shape[0],shape[1])


    else:
        dark_current = np.ones(shape) * ncc.cfg ['dark_current']

    ncc.create_var_auto(dims, dark_current, attr, 'f8')
