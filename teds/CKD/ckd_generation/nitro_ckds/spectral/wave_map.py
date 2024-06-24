"""
Wavelength Map with Spectral Smile
"""

import numpy as np
from scipy.interpolate import CubicSpline
from teds.CKD.ckd_generation.nitro_ckds.spectral.wavelength import import_spectral_smile_from_excel


def generate(ncc):
    row_ix, col_ix, act_2d, wl_2d = import_spectral_smile_from_excel(ncc.cfg)
    
    act = act_2d[:,0]
    wl = wl_2d[0,:]
    dims = ['detector_row', 'detector_col']
    nrows, ncols = ncc.get_shape(dims)
    rows = np.arange(nrows)
    cols = np.arange(ncols)
    
    # interpolate wavelengths over columns for each act pos
    # also interpolate row indices over columns
    wave_map1 = np.empty((len(act), ncols))
    row_ix1 = np.empty((len(act), ncols))
    for i in range(len(act)):
        spl_wl = CubicSpline(col_ix[i,:], wl)
        spl_row = CubicSpline(col_ix[i,:], row_ix[i,:])
        wave_map1[i] = spl_wl(cols)
        row_ix1[i] = spl_row(cols)
    
    # interpolate wavelengths over rows for each column
    wave_map2 = np.empty((nrows, ncols))
    for j in range(len(cols)):
        spl = CubicSpline(row_ix1[:,j],wave_map1[:,j])
        wave_map2[:,j] = spl(rows)
    
    attrs = {
        'long_name' : "wavelength map",
        'units' : 'nm'
    }
    ncc.create_var_auto(dims, wave_map2, attrs, 'f8')
