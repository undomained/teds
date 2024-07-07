"""
Wavelength
"""
import numpy as np

import pandas as pd
from scipy.interpolate import bisplrep, bisplev

def generate(ncc):
    cfg = ncc.cfg
    dims = ['across_track', 'detector_column']
    row_ix, col_ix, act_2d, wl_2d = import_spectral_smile_from_excel(cfg)
    wave_spline = bisplrep(act_2d, col_ix, wl_2d)
    act = act_2d[:,0]
    act_new = np.linspace(np.min(act), np.max(act), cfg['dimensions']['across_track'])
    cols = np.arange(0, cfg['dimensions']['detector_column'])
    wavelength = bisplev(act_new, cols, wave_spline)
    
    ncc.create_var_auto(dims, wavelength, dtype = 'f8')



def import_spectral_smile_from_excel(cfg):
    filename = cfg['paths']['dir_external'] + 'TANGO_Nitro_011_InstrumentModelInputs_20240307.xlsx'
    # Get centroid X position on detector (to be converted to columns)
    df_x = pd.read_excel(filename, 'SpectralSmile', skiprows = 6, nrows=21,  usecols= 'B:M')
    # Get centroid Y position on detector (to be converted to rows)
    df_y = pd.read_excel(filename, 'SpectralSmile', skiprows = 6, nrows = 21, usecols= 'O:U')
    
    # Trying to import horrible excel file...
    obsolete_columns = [col for col in df_x.columns if 'delta' in str(col)]
    obsolete_columns.append('Centroid')
    df_x.drop(obsolete_columns, axis = 'columns', inplace=True)
    df_y.columns = [colname.replace(".1","") for colname in df_y.columns]
    df_y.drop('Centroid', axis = 'columns', inplace=True)
    
    wavelengths =  [col for col in df_x.columns if 'ACT' not in str(col)]
    px_size = cfg['detector_pixel_size'] 
    act =  np.array(df_x['ACT Pos'])
    col_ix = np.array(df_x[np.array(wavelengths)])/px_size
    row_ix = np.array(df_y[np.array(wavelengths, dtype= str)])/px_size
    
    # assume data is mirrored with axis at nadir (actpos = 0)
    act = mirror_on_first_row(act, neg = True)
    row_ix = mirror_on_first_row(row_ix, neg = True)
    col_ix = mirror_on_first_row(col_ix)
    
    # assume image is projected on center of detector
    row_ix += cfg['dimensions']['detector_row']/2
    col_ix += cfg['dimensions']['detector_column']/2 
    
    act_2d = np.tile(act, (col_ix.shape[1],1)).T
    wl_2d  = np.tile(wavelengths, (col_ix.shape[0],1))
    
    return row_ix, col_ix, act_2d, wl_2d


def mirror_on_first_row(arr, neg = False):
    sign = -1 if neg else 1
    nrows = arr.shape[0]
    if len(arr.shape) > 1:
        ncols = arr.shape[1]
        mirrored = np.empty((2*nrows-1, ncols))
    else:
        mirrored = np.empty((2*nrows-1)) 
    mirrored[nrows-1:] = arr
    mirrored[:nrows-1] = sign * np.flip(arr[1:], axis = 0)
    return mirrored

