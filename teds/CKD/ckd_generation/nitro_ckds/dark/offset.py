"""
Dark Current
"""

import numpy as np


def generate(ncc):
    attr = {
        'long_name': 'detector dark offset',
        'units': 'counts',
        'comment': 'Placeholder, no data available'
    }
    
    dims = ['detector_row', 'detector_column']
    dark_offset = np.zeros(ncc.get_shape(dims))    
    ncc.create_var_auto(dims, dark_offset, attr, 'f8')
