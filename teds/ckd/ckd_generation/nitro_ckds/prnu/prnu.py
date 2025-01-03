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
    seed = cfg['seed']
    prnu_max = 1 + prnu_budget/100
    prnu_min = 1 - prnu_budget/100    

    QE = ncc.get_var('prnu/quantum_efficiency')

    attr = {
        'long_name': 'pixel response non-uniformity',
        'comment' :  'PRNU with Quantum Efficiency'
    }
    if seed is not None:
        np.random.seed(seed)
    prnu = np.random.uniform(prnu_min, prnu_max, size = ncc.get_shape(dims))

    prnu_QE = prnu * QE

    ncc.create_var_auto(dims, prnu_QE, attr, 'f8')

