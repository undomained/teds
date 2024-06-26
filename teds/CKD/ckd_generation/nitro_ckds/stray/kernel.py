"""
Straylight kernel variables
"""

import numpy as np
import h5py as h5
from scipy.interpolate import CubicSpline

def generate(ncc):
    cfg = ncc.cfg

    # Import data from carbon instrument
    with h5.File(f"{cfg['paths']['dir_external']}ckd_stray.nc") as f:
        kernel_rows = f['stray/kernel_rows'][:]
        kernel_cols = f['stray/kernel_cols'][:]
        kernel_fft_sizes = f['stray/kernel_fft_sizes'][:]
        kernels_fft = f['/stray/kernels_fft'][:]
        eta = f['stray/eta'][:]
        weights = f['stray/weights'][:]
        edges = f['stray/edges'][:]
    
    # reshape some variables for the time being
    C = cfg['dimensions']['detector_col']
    R = cfg['dimensions']['detector_row']
    col_new = np.arange(C)/C
    row_new = np.arange(R)/R
    col_old = np.arange(eta.shape[1])/eta.shape[1]
    row_old = np.arange(eta.shape[0])/eta.shape[0]
    eta_col_new = np.zeros((len(row_old), C))
    wei_col_new = np.zeros((len(weights), len(row_old), C))
    for r in range(len(row_old)):
        spl = CubicSpline(col_old, eta[r])
        eta_col_new[r] = spl(col_new)
        for i in range(len(weights)):
            splw = CubicSpline(col_old, weights[i,r,:])
            wei_col_new[i,r,:] = splw(col_new)

    eta_new = np.zeros((R, C))
    wei_new = np.zeros((len(weights), R, C))
    for c in range(len(col_new)):
        spl = CubicSpline(row_old, eta_col_new[:,c])
        eta_new[:,c] = spl(row_new)
        for i in range(len(weights)):
            splw = CubicSpline(row_old, wei_col_new[i,:,c])
            wei_new[i,:,c] = splw(row_new)

    # Create rows
    attr = {'long_name': "number of rows in each kernel"}
    ncc.create_var('kernel_rows', ['kernel'], kernel_rows, attr, 'u1')
    # Create cols
    attr = {'long_name': "number of columns in each kernel"}
    ncc.create_var('kernel_cols', ['kernel'], kernel_cols, attr, 'u1')
    # Create fft_sizess
    attr = {'long_name': "sizes of kernel FFTs"}
    ncc.create_var('kernel_fft_sizes', ['kernel'], kernel_fft_sizes, attr, 'u1')
    # Create ffts
    attr = {'long_name': 'Fourier transforms of the kernels', 'units':'1'}
    ncc.create_var('kernels_fft', ['kernel_fft_size'], kernels_fft, attr, 'f8')
    # create eta
    dims_img = ['detector_row', 'detector_col']
    eta_attr = {'long_name': 'internal scattering factor'}
    ncc.create_var('eta', dims_img, eta_new, eta_attr, 'f8')
    # create weights
    dims_weights = np.append('kernel', dims_img)
    attr = {'long_name': 'kernel weights', 'units':'1'}
    ncc.create_var('weights', dims_weights, wei_new, attr, 'f8')
    # create edges
    dims_edges = ['kernel', 'edges_of_box']
    attr = {'long_name': 'distances of subimage edges from the detector edges'}
    ncc.create_var('edges', dims_edges, edges, attr, 'f8')
    









    
