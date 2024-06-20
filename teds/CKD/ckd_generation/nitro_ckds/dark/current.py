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
    
    dims = ['detector_row', 'detector_col']
    dark_current = np.random.uniform(10, 360, size = ncc.get_shape(dims))
    ncc.create_var_auto(dims, dark_current, attr)
