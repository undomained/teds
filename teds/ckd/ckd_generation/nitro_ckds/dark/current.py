"""
Dark Current
"""

import numpy as np

def generate(ncc):
    attr = {
        'long_name': 'detector dark current',
        'units': 'counts/s',
        'comment': 'Not temperature dependent'
    }
    
    dims = ['detector_row', 'detector_column']
    shape = ncc.get_shape(dims)
    dark_current = np.ones(shape) * ncc.cfg ['dark_current']
    
    ncc.create_var_auto(dims, dark_current, attr, 'f8')
