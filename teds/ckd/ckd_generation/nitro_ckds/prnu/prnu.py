"""
Pixel Response Non-Uniformity
Does not include Quantum Efficiency
"""

import numpy as np

from scipy.ndimage import gaussian_filter

def generate(ncc):
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']
    prnu_budget = cfg['prnu'] # in %, from TANGO Nitro SNR budget (TNO)
    prnu_max = 1 + prnu_budget/100
    prnu_min = 1 - prnu_budget/100    

    attr = {
        'long_name': 'pixel response non-uniformity',
        'comment' :  'Placeholder, does not include Quantum Efficiency'
    }
    prnu = np.random.uniform(prnu_min, prnu_max, size = ncc.get_shape(dims))
    ncc.create_var_auto(dims, prnu, attr, 'f8')

