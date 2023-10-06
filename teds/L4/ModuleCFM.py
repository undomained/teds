#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 2023.

@author: Manu Goudar
"""

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator as RGI
from shapely.geometry import LineString
from scipy.ndimage import gaussian_filter1d

from .ModuleDataContainer import DataCont
from .ModulePlume import get_thresholdbasedplume
from .ModuleRemoveBackground import remove_background


def rearrange_data(x, y, _ix, _nx):
    # these 0.35 and 3 are based on the pixels before the
    # source pixel or after the source pixels
    # sometimes we have plume pixels (before the plume too)
    if (((_nx - _ix) / _nx) < 0.35) or ((_nx - _ix) < 3):
        x1 = x[:_ix+1]
        y1 = y[:_ix+1]
        # reset the order start with source
        return x1[::-1], y1[::-1]
    else:
        # reset the order start with source
        return x[_ix-1:], y[_ix-1:]


def fit_plumeline(x, y, sx):
    _nx = len(x)
    # fit polynominal
    _ix = np.searchsorted(x, sx)
    x1, y1 = rearrange_data(x, y, _ix, _nx)
    sgm = np.minimum(len(x1)*0.05, 5)
    return x1, gaussian_filter1d(y1, sigma=sgm, mode='nearest')


def getplumeline(df, src):
    """Extract a plume line

    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe containing x,y,conc
    src : Vector{Float64}
        Source location
    """
    if df.x.std() > df.y.std():
        di = (df.sort_values(by='j')).groupby(by="j").apply(
            lambda f: [np.mean(f['x']), sum(f['y']*f['conc'])/f['conc'].sum()])
        dfi = pd.DataFrame(di.to_list(), columns=['x', 'y'])
        xvi, yvi = fit_plumeline(dfi.x.values, dfi.y.values, src[0])
        return np.vstack([dfi.x.values, dfi.y.values]).T, np.vstack([xvi, yvi]).T
    else:
        di = (df.sort_values(by='i')).groupby(by="i").apply(
            lambda f: [np.mean(f['y']), sum(f['x']*f['conc'])/f['conc'].sum()])
        dfi = pd.DataFrame(di.to_list(), columns=['y', 'x'])
        yvi, xvi = fit_plumeline(dfi.y.values, dfi.x.values, src[1])
        return np.vstack([dfi.x.values, dfi.y.values]).T, np.vstack([xvi, yvi]).T


"""
functions to segment the plume line and draw transaction lines
"""


def segment_plumeline(segment_ds, fitted_pline):
    """Segment line at equal distances

    Segments a line into equal distances based on the sample distance

    Parameters
    ----------
    segment_ds : Float
        Sampling distance
    fitted_pline : Array (n,2)
        Coordinates representing a line

    Examples
    --------
    FIXME: Add docs.

    """
    # Convert given coordinates to a Shapely linestring
    fpline = LineString(tuple(map(tuple, fitted_pline)))
    # compute number of points based on the length of the line
    n = np.int_(fpline.length/segment_ds)
    # segment the line based on sampling distance
    segmented_line = LineString([(fpline.interpolate(i * segment_ds)) for i in range(n)])
    # Convert the linestring to coordinates array
    plumeline = (np.array(segmented_line.xy)).T
    # compute direction vector. Direction at each point
    vec = np.zeros_like(plumeline)
    vec[:-1] = (plumeline[1:] - plumeline[:-1])
    vec[-1, :] = vec[-2, :]
    direction_vector = (vec.T/np.linalg.norm(vec, axis=1)).T
    return n, segmented_line, direction_vector


def get_transactionlines(n, ds, length_side, pline):
    """Compute transaction lines

    A set parallel lines to the plumeline are dran and sampled. Then these points are rearranged to
    get the coordinates of transaction lines.

    Parameters
    ----------
    n : Int
        Number of transaction lines
    ds : Float
        Sampling distance in a transaction line
    length_side : Float
        Lenth of a side of transaction line
    pline : Shapely Linestring
        Plumeline that is segmented as 'n' points

    Examples
    --------
    FIXME: Add docs.

    """
    # total number of points on transaction line
    m = 2*np.int_(length_side//ds) + 1
    # number of points on one side of the transaction line
    m1 = np.int_((m - 1)/2)
    # total number of points on the transaction line
    line_pts = np.zeros((n, m, 2))
    for i in range(m1):
        left = pline.parallel_offset(ds*(i+1), 'left')
        right = pline.parallel_offset(ds*(i+1), 'right')
        line_pts[:, m1, :] = np.array(pline.coords)
        dxl = left.length/(n-1)
        dxr = right.length/(n-1)
        for j in range(n):
            line_pts[j, m1-1-i, :] = np.array(left.interpolate((j) * dxl).coords)[0]
            line_pts[j, m1+1+i, :] = np.array(right.interpolate((j) * dxr).coords)[0]
    return line_pts, np.arange(-length_side, length_side+1, ds)


def get_tlines(fittedline, pline_ds, tline_ds, tline_sidelen, interp):
    """Segment plume line and compute transaction line

    Parameters
    ----------
    fittedline : Vector (n,2)
        A Line segment representing plumeline
    pline_ds : Float64
        Sampling distance between transaction lines
    tline_ds : Float
        Sampling distance in a transaction line
    side_length : Float
        Distance of the transaction line on one side of the plume line
    interp : Interpolation function
        Interpolation function to interpolate concentration values

    Examples
    --------
    FIXME: Add docs.

    """
    # segment the fitted plume lines based on a distance 'segmentplumeline_ds'
    n, plumeline, direction_vec = segment_plumeline(pline_ds, fittedline)
    pline_coords = np.array(plumeline.coords.xy).T

    # Transactional lines.
    # Draw parallel lines to the plumeline and segment it
    _lines, xdata = get_transactionlines(n, tline_ds, tline_sidelen, plumeline)

    tlines = []
    for i in range(n):
        _ln = DataCont()
        _ln.__setattr__('directional_vector', direction_vec[i])
        _ln.__setattr__('dist_from_src', (i+1)*pline_ds)

        # initial data
        i_data = DataCont()
        i_data.__setattr__('coords', _lines[i, :, :])
        ydata = interp(i_data.coords)
        i_data.__setattr__('co', ydata)
        i_data.__setattr__('line', xdata)
        _ln.__setattr__("initialdata", i_data)      # add data to line
        # remove background
        backdata = remove_background(_ln.initialdata)
        # define final variables if remove_background is successful
        if backdata.backgroundremovalflag:
            _ix = backdata.filtered_index
            _ln.__setattr__("co", backdata.back_removed[_ix])
            _ln.__setattr__("line", backdata.cutoff_line[_ix])
            _ln.__setattr__("coords", backdata.cutoff_coords[_ix, :])
        tlines.append(_ln)
    return pline_coords, tlines


def get_emission(tlines, interp_u, interp_v, dxy):
    """Estimate emission

    Parameters
    ----------
    tlines : List
        Transaction lines
    interp_u : Interpolation
        Interpolation function for u velocity
    interp_v : Interpolation
        Interpolation function for v velocity
    dxy : Float
        Sampling distance of the transaction line

    """
    emis = []
    for _ln in tlines:
        vel_mag = (interp_u(_ln.coords)*_ln.directional_vector[0]
                   + interp_v(_ln.coords)*_ln.directional_vector[1])
        _ln.__setattr__("vel_mag", vel_mag)
        emission = np.sum(_ln.vel_mag*_ln.co)*dxy
        _ln.__setattr__('emission', emission)
        emis.append(emission)
    return emis


def get_massflux(conc, grid, cdfthreshweight, interp_u, interp_v,
                 dist_tline=100, tline_ds=50, tline_sidelen=5000):
    """Compute cross-sectional mass flux

    Parameters
    ----------
    conc : Matrix
        Concentration of a gas in kg/m^2
    grid : Grid class
        Class describing a grid structure
    cdfthreshweight : Float
        Avalue of cumulative distribution function which is used to compute threshold
    interp_u : Interpolation
        Interpolation function for u velocity
    interp_v : Interpolation
        Interpolation function for v velocity
    dist_tline : Float, optional
        Distance between Transaction lines. Defaults to 100m
    tline_ds : Float, optional
        Sampling distance of a transaction line. Defaults to 50m
    tline_sidelen : Float, optional
        Length of a side of a transaction line. Defaults to 5000m

    Returns
    --------
    massflux : class
        A class containing all data for emission estimation
    emission : Float
        Estimated emission

    Examples
    --------
    FIXME: Add docs.

    """

    # find the pixel with the origin
    x = np.sqrt((grid.xc - grid.sourcexy[0])**2 +
                (grid.yc - grid.sourcexy[1])**2)
    origin_id = np.unravel_index(np.argmin(x), x.shape)

    # detect plume
    massflux = get_thresholdbasedplume(conc, origin_id, cdfthreshweight)

    if massflux.flag:
        pts = np.column_stack((grid.xc.ravel(), grid.yc.ravel()))
        # create interpolation for data
        interp = RGI(pts, conc.ravel(), fill_value=0)

        # compute plumeline and smoothen it
        i, j = np.where(massflux.plume)
        df_pline = pd.DataFrame.from_dict({"i": i, "j": j,
                                           "x": grid.xc[massflux.plume],
                                           "y": grid.yc[massflux.plume],
                                           "conc": conc[massflux.plume]})
        rawpline, pline = getplumeline(df_pline, grid.sourcexy)

        # Segment plume line and compute transaction lines
        plumeline, tlines = get_tlines(pline, dist_tline, tline_ds,
                                       tline_sidelen, interp)

        # write data to massflux calss
        massflux.__setattr__("fitted_plumeline", pline)
        massflux.__setattr__("distancebetweentlines", dist_tline)
        massflux.__setattr__("tline_ds", tline_ds)
        massflux.__setattr__("segmentedplumeline", plumeline)
        massflux.__setattr__("tlines", tlines)

        emis = get_emission(massflux.tlines, interp_u,
                            interp_v, massflux.tline_ds)
        return massflux, np.mean(np.array(emis))
    else:
        return massflux, None
