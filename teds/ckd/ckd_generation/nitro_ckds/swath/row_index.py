"""
Row index
"""
import numpy as np

from scipy.interpolate import CubicSpline
from teds.ckd.ckd_generation.nitro_ckds.spectral.wavelength import import_spectral_smile_from_excel

def generate(ncc):
    cfg = ncc.cfg
    dims = ['across_track', 'detector_column']

    # import from external data excel
    row_ix, col_ix, act_2d, _ = import_spectral_smile_from_excel(cfg)
    act = act_2d[:,0]
    cols = np.arange(0, cfg['dimensions']['detector_column'])
    act_new = np.linspace(np.min(act), np.max(act), cfg['dimensions']['across_track'])

    # Since x and y is non-monotonic, interpolate each dimension seperately
    row_ix_intermediate = np.zeros((len(act), len(cols)))
    row_indices = np.zeros((len(act_new), len(cols)))
    for i, ap in enumerate(act):
        spl = CubicSpline(col_ix[i, :], row_ix[i, :])
        row_ix_intermediate[i,:] = spl(cols)
    for i, col in enumerate(cols):
        spl = CubicSpline(act, row_ix_intermediate[:, i])
        row_indices[:, i] = spl(act_new)

    # make a 2d spline 
    # los_spline = bisplrep(act_2d, col_ix, row_ix, kx = 3, ky = 3)
    # act_new = np.linspace(np.min(act), np.max(act), cfg['dimensions']['across_track'])
    # cols = np.arange(0, cfg['dimensions']['detector_column'])
    # row_indices = bisplev(act_new, cols, los_spline)

    
    attr = {
        'long_name': 'Row index at which to extract the spectra',
    }
    
    ncc.create_var_auto(dims, row_indices, attr, 'f8')
