"""
Straylight kernel variables
"""

import numpy as np
import h5py as h5
from scipy.interpolate import CubicSpline
#import urllib.request
import os
from getpass import getpass
import requests

    
def generate(ncc):
    cfg = ncc.cfg
    # Transfer data form ckd_stray to ckd with all data

    # Import data
    straydatapath = f"{cfg['paths']['dir_external']}ckd_stray_nitro.nc"
    straydatalink = "https://surfdrive.surf.nl/files/index.php/s/xhMKNPX5SoBAfA3/download?path=%2Fdata&files=ckd_stray_nitro.nc"
    if not os.path.isfile(straydatapath):
        print("{:=^60}".format(" Straylight CKD Data missing "))
        print("Download CKD straylight nitro files from TANGO surfdrive")
        print("url: {}".format(straydatalink))
        print("target_folder: {}".format(straydatapath))
        print("{:=^60}".format(""))
        exit()
    
    with h5.File(straydatapath) as f:
        kernel_rows = f['kernel_rows'][:]
        kernel_cols = f['kernel_cols'][:]
        kernel_fft_sizes = f['kernel_fft_sizes'][:]
        kernel_fft_re = f['kernel_fft_re'][:]
        kernel_fft_im = f['kernel_fft_im'][:]
        eta = f['eta'][:]
        weights = f['weights'][:]
        edges = f['edges'][:]

    # create dimensions for kernel fft size and number of kernels
    ncc.create_dim('coarse_image_rows', eta.shape[0])
    ncc.create_dim('coarse_image_cols', eta.shape[1])
    ncc.create_dim('kernel_fft_size', len(kernel_fft_re))
    ncc.create_dim('kernel', len(kernel_rows))
    ncc.create_dim('edges_of_box', 4)
    ncc.create_dim('n_van_cittert', cfg["n_van_cittert"])

    # Create rows
    attr = {'long_name': "number of rows in each kernel"}
    ncc.create_var('kernel_rows', ['kernel'], kernel_rows, attr, 'i4')
    # Create cols
    attr = {'long_name': "number of columns in each kernel"}
    ncc.create_var('kernel_cols', ['kernel'], kernel_cols, attr, 'i4')
    # Create fft_sizess
    attr = {'long_name': "sizes of kernel FFTs"}
    ncc.create_var('kernel_fft_sizes', ['kernel'], kernel_fft_sizes, attr, 'i4')
    # Create ffts
    attr = {'long_name': 'Real part of fourier transforms of the kernels', 'units':'1'}
    ncc.create_var('kernel_fft_re', ['kernel_fft_size'], kernel_fft_re, attr, 'f8')
    # Create ffts
    attr = {'long_name': 'Imaginary part of fourier transforms of the kernels', 'units':'1'}
    ncc.create_var('kernel_fft_im', ['kernel_fft_size'], kernel_fft_im, attr, 'f8')
    # create eta
    dims_img = ['coarse_image_rows', 'coarse_image_cols']
    eta_attr = {'long_name': 'internal scattering factor'}
    ncc.create_var('eta', dims_img, eta, eta_attr, 'f8')
    # create weights
    dims_weights = np.append('kernel', dims_img)
    attr = {'long_name': 'kernel weights', 'units':'1'}
    ncc.create_var('weights', dims_weights, weights, attr, 'f8')
    # create edges
    dims_edges = ['kernel', 'edges_of_box']
    attr = {'long_name': 'distances of subimage edges from the detector edges'}
    ncc.create_var('edges', dims_edges, edges, attr, 'f8')

    

    #ncc.nc.createVariable(|"n_van_cittert", "i4")



# def generate(ncc):
#     cfg = ncc.cfg
#     # Transfer data form ckd_stray to ckd with all data
#     straydatapath = f"{cfg['paths']['dir_external']}ckd_stray_nitro.nc"
    
#     with h5.File(straydatapath) as f:
#         kernel = f['kernel'][:]
#         crow = f['crow'][:]
#         ccol = f['ccol'][:]

#     print(kernel.shape)
#     # create dimensions for kernel fft size and number of kernels
#     ncc.create_dim('coarse_image_rows', kernel.shape[1])
#     ncc.create_dim('coarse_image_cols', kernel.shape[2])
#     ncc.create_dim('kernel', len(kernel))

#     # Create rows
#     attr = {'long_name': "row index of kernel center"}
#     ncc.create_var('crow', ['kernel'], crow, attr, 'i4')
#     # Create cols
#     attr = {'long_name': "column index of kernel center"}
#     ncc.create_var('ccol', ['kernel'], ccol, attr, 'i4')
#     # Create kernels
#     attr = {'long_name': "kernels"}
#     ncc.create_var('kernel', ['kernel', 'coarse_image_rows', 'coarse_image_cols'], kernel, attr, 'f8')
    









    
