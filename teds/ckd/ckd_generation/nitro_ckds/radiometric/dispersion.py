import pandas as pd
from scipy.interpolate import CubicSpline
import numpy as np

def generate(ncc):
    """
    Dispersion
    """
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']
    attr = {
        'long_name': 'Dispersion',
        'comment' : '',
        'units' : 'um'
    }

    # create dispersion variable
    dispersion = import_dispersion_from_excel(ncc, cfg)
    ncc.create_var_auto(dims, dispersion, attr, 'f8')

def import_dispersion_from_excel(ncc, cfg):
    filename = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_InstrumentModelInputs_20240307.xlsx'
    # Get dispersion
    df = pd.read_excel(filename, 'Dispersion', skiprows = 7, nrows=7,  usecols= 'B:G', index_col = 0)
    wl_in = np.array(df.columns)
    ap_in = np.array(df.index)
    ap_out = np.linspace(ap_in[0], ap_in[-1], cfg["dimensions"]["across_track"])

    # First interpolate ACT input dimension 
    disp_to_ap = np.zeros((cfg["dimensions"]["across_track"], len(wl_in)), dtype = float) # (ACT, 5)
    for j, wl in enumerate(wl_in):
        disp = df.iloc[:, j] # dispersion per wavelength over input ACT
        spl = CubicSpline(ap_in, disp)
        disp_to_ap[:, j] = spl(ap_out)

    # Get wavelengths per ACT, and interpolate dispersion to those wavelength
    wavelength_per_actpos = ncc.get_var('spectral/wavelength')
    disp_to_col = np.zeros((cfg["dimensions"]["across_track"], cfg["dimensions"]["detector_column"]), dtype = float) # (ACT, col)
    for i, ap in enumerate(ap_out):
        disp = disp_to_ap[i, :]
        spl = CubicSpline(wl_in, disp)
        wl_out = wavelength_per_actpos[i, :]
        disp_to_col[i, :] = spl(wl_out)

    # Project ACT to rows
    row_indices = ncc.get_var('swath/row_index')
    row_out = np.arange(cfg["dimensions"]["detector_row"])
    disp_to_det = np.zeros((cfg["dimensions"]["detector_row"], cfg["dimensions"]["detector_column"]), dtype = float)
    for j in range(cfg["dimensions"]["detector_column"]):
        disp = disp_to_col[:, j]
        row_this_col = row_indices[:, j]
        spl = CubicSpline(row_this_col, disp)
        disp_to_det[:, j] = spl(row_out)

    return disp_to_det

