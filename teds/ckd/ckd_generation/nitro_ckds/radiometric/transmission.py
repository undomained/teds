import pandas as pd
from scipy.interpolate import CubicSpline
import numpy as np

def generate(ncc):
    """
    transmission
    """
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']
    attr = {
        'long_name': 'Transmission',
        'comment' : 'fraction',
        'units' : '1'
    }

    # create transmission variable
    transmission = import_transmission_from_excel(ncc, cfg)
    ncc.create_var_auto(dims, transmission, attr, 'f8')

def import_transmission_from_excel(ncc, cfg):
    filename = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_InstrumentModelInputs_20240307.xlsx'
    # Get transmission from M00 term of Mueller matrix
    df = pd.read_excel(filename, 'MuellerMat_ManufGratingFinal', skiprows = 5, nrows=127,  usecols= 'A:B,H,N', index_col = 0)
    df.columns = ["act_pos 0.0", "act_pos 0.7", "act_pos 1.0"]
    ap_in = [0, 0.7, 1.0]

    wl_in = np.array(df.index)
    ap_out = np.linspace(ap_in[0], ap_in[-1], cfg["dimensions"]["across_track"])

    # First interpolate ACT input dimension 
    tr_to_ap = np.zeros((cfg["dimensions"]["across_track"], len(wl_in)), dtype = float) # (ACT, wl_in)
    for j, wl in enumerate(wl_in):
        tr = df.iloc[j, :] # transmission per wavelength over input ACT
        spl = CubicSpline(ap_in, tr)
        tr_to_ap[:, j] = spl(ap_out)

    # Get wavelengths per ACT, and interpolate transmission to those wavelength
    wavelength_per_actpos = ncc.get_var('spectral/wavelength')
    tr_to_col = np.zeros((cfg["dimensions"]["across_track"], cfg["dimensions"]["detector_column"]), dtype = float) # (ACT, col)
    for i, ap in enumerate(ap_out):
        tr = tr_to_ap[i, :]
        spl = CubicSpline(wl_in, tr)
        wl_out = wavelength_per_actpos[i, :]
        tr_to_col[i, :] = spl(wl_out)

    # Project ACT to rows
    row_indices = ncc.get_var('swath/row_index')
    row_out = np.arange(cfg["dimensions"]["detector_row"])
    tr_to_det = np.zeros((cfg["dimensions"]["detector_row"], cfg["dimensions"]["detector_column"]), dtype = float)
    for j in range(cfg["dimensions"]["detector_column"]):
        tr = tr_to_col[:, j]
        row_this_col = row_indices[:, j]
        spl = CubicSpline(row_this_col, tr)
        tr_to_det[:, j] = spl(row_out)

    return tr_to_det

