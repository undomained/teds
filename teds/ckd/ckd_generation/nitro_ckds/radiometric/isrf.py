"""
Instrument Spectral Response Function
"""
import numpy as np

import os
import csv
from scipy.interpolate import CubicSpline
from teds.ckd.ckd_generation.nitro_ckds.spectral.wavelength import mirror_on_first_row

def generate(ncc): # Get the data and combine
    cfg = ncc.cfg
    data_dir = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_ISRF_SonyTolerancedMp2S_LowResSampling/'
    data_files = os.listdir(data_dir)

    wl, ap = [], []  #wavelengths, act_pos
    # Get files, as well as wavelength and actpos data (from filenames...)
    for fi, f in enumerate(data_files):
        parts = f.split('_')
        wl.append(float(parts[-2][:-2])) # wavelength in nm
        ap.append(int(parts[-1][:-10])) # across track position 

    # Combine data in dataset
    wl_in, ap_in = np.unique(wl), np.unique(ap)
    for i, file in enumerate(data_files):
        with open(data_dir + file, 'r', newline='\n') as f:
            rows = csv.reader(f, delimiter = ',')
            filedata = np.array([row for row in rows])
            filedata = filedata[1:, :]
            
            if i == 0:  # initialize data array
                data = np.zeros((len(wl_in), len(ap_in), len(filedata), len(filedata[0])))
            p = np.flatnonzero(wl_in == wl[i])
            q = np.flatnonzero(ap_in == ap[i])
            data[p,q,:,:] = filedata
            # p = wavelength
            # q = across track position
            # gives a matrix containing columns:
            # physical distance um, spectral distance nm, ISRF
    
    # Interpolate and Calculate ISRF
    isrf_in = data[:,:,:,2] # isrf data
    d_spectral_in = data[:,:,:,1] # spectral distance of data
    d_physical_in = data[:,:,:,0] # physical distance of data

    min_wl = cfg['lbl_min_wl']
    max_wl = cfg['lbl_max_wl']
    n_wl = cfg['dimensions']['lbl_samples']
    res = (max_wl - min_wl) / (n_wl - 1)

    d_spectral_out = np.arange(0, 0.8 + res, step = res)
    d_spectral_out = mirror_on_first_row(d_spectral_out, neg = True)

    X = n_wl  # samples in spectral direction
    Y = cfg['dimensions']['across_track'] # samples in spatial direction
    N = len(d_spectral_out)  # output number of samples of ISRF shape

    ap_out = np.linspace(-1, 1, Y)
    # 1. convert all isrf's to the same x-axis (spectral distance)
    isrf_1 = np.empty((len(wl_in), len(ap_in), N))
    isrf_2 = np.empty((len(wl_in), Y, N))
    for i in range(len(wl_in)):
        for j in range(len(ap_in)):
            spl = CubicSpline(d_spectral_in[i,j], isrf_in[i,j])
            isrf_1[i, j] = spl(d_spectral_out)
        # 2. interpolate isrf's over across track positions
        for k in range(len(d_spectral_out)):
            # only two datapoints so cubicspline might be overkill
            spl = CubicSpline(ap_in, isrf_1[i,:,k])
            isrf_2[i,:,k] = spl(ap_out)

    # 3. interpolate isrf's over wavelengths
    wl_out = np.linspace(cfg['lbl_min_wl'], cfg['lbl_max_wl'], X)
    isrf_out = np.empty((X,Y,N))
    for j in range(len(ap_out)):
        for k in range(len(d_spectral_out)):
            spl = CubicSpline(wl_in, isrf_2[:,j,k])
            isrf_out[:,j,k] = spl(wl_out) 
    
    # Array with wavelength range per isrf
    isrf_wl = np.empty((X,Y,N)) 
    for i, wl in enumerate(wl_out):
        for j in range(len(ap_out)):
            isrf_wl[i, j] = d_spectral_out + wl
    
    isrf_out = np.transpose(isrf_out,(1, 0, 2))
    isrf_wl = np.transpose(isrf_wl, (1, 0, 2))

    # Create ISRF ckd
    newdims = {'isrf_samples' : N}
    ncc.create_dims_auto(newdims)
    dims = ['across_track', 'lbl_samples', 'isrf_samples']
    attr = {
        'source': 'TANGO_Nitro_011_ISRF_SonyTolerancedMp2S_LowResSampling (TNO)',
        'comment' : '[!] Make sure the isrf_sample increments agree with the increments in the line-by-line spectra for correction convolution'
    }
    ncc.create_var_auto(dims, isrf_out, attr, dtype = 'float32')

    attr_wl = {
        'source': 'TANGO_Nitro_011_ISRF_SonyTolerancedMp2S_LowResSampling (TNO)',
        'comment': 'wavelength values per isrf', 
        'units' : 'nm'
    }
    ncc.create_var('isrf_wavelengths', dims, isrf_wl, attr_wl, dtype = 'float32')
