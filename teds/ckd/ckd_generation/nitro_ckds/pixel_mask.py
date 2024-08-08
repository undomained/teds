"""
Detector Pixel Mask
"""
import numpy as np

def generate(ncc):
    dims = ['detector_row', 'detector_column']
    mask = np.zeros(ncc.get_shape(dims))

    attr = {
        'long_name': 'detector pixel mask',
        'valid_range' : (np.int8(0), np.int8(1)),
        'flag_values' : (np.int8(0), np.int8(1)),
        'flag_meanings' : '0:good, 1:bad'
        }

    ncc.create_var_auto(dims, mask, attr)
