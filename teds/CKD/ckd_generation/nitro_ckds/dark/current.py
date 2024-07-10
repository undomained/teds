"""
Dark Current
"""

import numpy as np

def generate(ncc):
    attr = {
        'long_name': 'detector dark current',
        'units': 'counts/s',
        'comment': 'Placeholder, random noise between BOL and EOL values'
    }
    
    dims = ['detector_row', 'detector_column']
    np.random.seed(100)
#    dark_current = np.random.uniform(10, 360, size = ncc.get_shape(dims))
    # Try something with a tail
#    mean_dc = 300.
#    sigma_dc = 100.
    mean_dc = 50.
    sigma_dc = 50.
    normal_std = np.sqrt(np.log(1 + (sigma_dc/mean_dc)**2))
    normal_mean = np.log(mean_dc) - normal_std**2 / 2
    dark_current =  np.random.lognormal(normal_mean, normal_std, size=ncc.get_shape(dims))
    ncc.create_var_auto(dims, dark_current, attr, 'f8')
