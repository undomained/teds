"""
Spatial Sample to Pixel Conversion
"""

import numpy as np


def generate(ncc):
    attr =  {
        'long_name': 'line of sight vectors' ,
        'comment' : 'contains dummy vector, will be used for geolocation'
    }
    dims = ['across_track',  'vector']
    los_vectors = np.tile([1, 0, 0], (ncc.get_shape(dims)[0], 1) )
    ncc.create_var_auto(dims, los_vectors, attr)
