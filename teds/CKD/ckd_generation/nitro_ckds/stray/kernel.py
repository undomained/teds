"""
Straylight kernel variables
"""

import numpy as np

def generate(ncc):
    cfg = ncc.cfg

    # create kernel_rows
    # do we need this?
    
    # create kernels_cols
    # do we need this?

    # create kernel_fft_sizes
    
    # create eta
    dims_img = ['detector_row', 'detector_col']
    eta = np.ones(ncc.get_shape(dims_img))
    eta_attr = {'long_name': 'internal scattering factor'}
    ncc.create_var('eta', dims_img, eta, eta_attr, 'f8')

    # create weights
    dims_weights = np.append('kernel', dims_img)
    weights = np.ones(ncc.get_shape(dims_weights))
    weights_attr = {'long_name': 'kernel weights'}
    ncc.create_var('weights', dims_weights, weights, weights_attr, 'f8')

    # create edges
    dims_edges = ['kernel', 'edges_of_box']
    edges = np.ones(ncc.get_shape(dims_edges))
    edges_attr = {'long_name': 'distances of subimage edges from the detector edges'}
    ncc.create_var('edges', dims_edges, edges, edges_attr, 'f8')




