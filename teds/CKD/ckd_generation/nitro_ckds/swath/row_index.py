"""
Row index
"""
import numpy as np

from scipy.interpolate import bisplrep, bisplev
from teds.CKD.ckd_generation.nitro_ckds.spectral.wavelength import import_spectral_smile_from_excel

def generate(ncc):
    cfg = ncc.cfg
    dims = ['across_track', 'detector_col']

    # import from external data excel
    row_ix, col_ix, act_2d, _ = import_spectral_smile_from_excel(cfg)
    act = act_2d[:,0]

    # make a 2d spline 
    los_spline = bisplrep(act_2d, col_ix, row_ix, kx = 3, ky = 3)
    act_new = np.linspace(np.min(act), np.max(act), cfg['dimensions']['across_track'])
    cols = np.arange(0, cfg['dimensions']['detector_col'])
    row_indices = bisplev(act_new, cols, los_spline)

    
    attr = {
        'long_name': 'Row index at which to extract the spectra',
        'assumption1': 'data is mirrored at nadir axis',
        'assumption2': 'image is projected on center of detector'
    }
    
    ncc.create_var_auto(dims, row_indices, attr, 'f8')
