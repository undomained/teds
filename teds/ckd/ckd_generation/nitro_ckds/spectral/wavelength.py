"""
Wavelength
"""
import numpy as np

import pandas as pd
from scipy.interpolate import bisplrep, bisplev
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import CubicSpline


def generate(ncc):
    cfg = ncc.cfg
    dims = ['across_track', 'detector_column']
    row_ix, col_ix, act_2d, wl_2d = import_spectral_smile_from_excel(cfg)
    act = act_2d[:,0]
    act_new = np.linspace(np.min(act), np.max(act), cfg['dimensions']['across_track'])
    cols = np.arange(0, cfg['dimensions']['detector_column'])

    # Since x and y is non-monotonic, interpolate each dimension seperately
    wavelength_intermediate = np.zeros((len(act), len(cols)))
    wavelength = np.zeros((len(act_new), len(cols)))
    for i, ap in enumerate(act):
        spl = CubicSpline(col_ix[i, :], wl_2d[i, :])
        wavelength_intermediate[i,:] = spl(cols)
    for i, col in enumerate(cols):
        spl = CubicSpline(act, wavelength_intermediate[:, i])
        wavelength[:, i] = spl(act_new)
    
    ncc.create_var_auto(dims, wavelength, dtype = 'f8')

def import_spectral_smile_from_excel(cfg):
    filename = cfg['paths']['im_inputs_tno']
    # Get centroid positions, x = act direction, y = spectral direction
    df = pd.read_excel(filename, "Centroids", skiprows = 6, nrows = 541, usecols= "B:F")
    df.columns = ['wavelength', 'act', 'alt', 'x', 'y']
    df['wavelength'] *= 1000 # convert wavelengths to nm
    wl =  np.unique(df['wavelength']) 
    act = np.unique(df['act'])
    alt = np.unique(df['alt'])

    # convert centroid position to pixel index
    center_pixel_row = cfg['dimensions']['detector_row']/2
    center_pixel_col = cfg['dimensions']['detector_column']/2
    df['x'] = df['x'] / cfg['detector_pixel_size'] + center_pixel_row
    df['y'] = df['y'] / cfg['detector_pixel_size'] + center_pixel_col

    row_ix = np.zeros((len(act), len(wl)))
    col_ix = np.zeros((len(act), len(wl)))
    act_2d = np.zeros((len(act), len(wl)))
    wl_2d  = np.zeros((len(act), len(wl)))
    for i, ap in enumerate(act):
        df_ap = df[(df['act'] == ap) & (df['alt'] == 0)] # using only center ALT position in slit
        row_ix[i,:] = df_ap["x"]
        col_ix[i,:] = df_ap["y"]
        act_2d[i,:] = df_ap["act"]
        wl_2d[i,:] = df_ap["wavelength"]

    # TODO: wachten op antwoord TNO over de illuminated area dimensions. 

    # df_x = pd.read_excel(filename, 'SpectralSmile', skiprows = 6, nrows=21,  usecols= 'B:M')
    # # Get centroid Y position on detector (to be converted to rows)
    # df_y = pd.read_excel(filename, 'SpectralSmile', skiprows = 6, nrows = 21, usecols= 'O:U')
    
    # # Trying to import horrible excel file...
    # obsolete_columns = [col for col in df_x.columns if 'delta' in str(col)]
    # obsolete_columns.append('Centroid')
    # df_x.drop(obsolete_columns, axis = 'columns', inplace=True)
    # df_y.columns = [colname.replace(".1","") for colname in df_y.columns]
    # df_y.drop('Centroid', axis = 'columns', inplace=True)
    
    # wavelengths =  [col for col in df_x.columns if 'ACT' not in str(col)]
    # px_size = cfg['detector_pixel_size'] 
    # act =  np.array(df_x['ACT Pos'])
    # col_ix = np.array(df_x[np.array(wavelengths)])/px_size
    # row_ix = np.array(df_y[np.array(wavelengths, dtype= str)])/px_size
    
    # # assume data is mirrored with axis at nadir (actpos = 0)
    # act = mirror_on_first_row(act, neg = True)
    # row_ix = mirror_on_first_row(row_ix, neg = True)
    # col_ix = mirror_on_first_row(col_ix)
    
    # # assume image is projected on center of detector
    # row_ix += cfg['dimensions']['detector_row']/2
    # col_ix += cfg['dimensions']['detector_column']/2 
    
    # act_2d = np.tile(act, (col_ix.shape[1],1)).T
    # wl_2d  = np.tile(wavelengths, (col_ix.shape[0],1))
    
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

