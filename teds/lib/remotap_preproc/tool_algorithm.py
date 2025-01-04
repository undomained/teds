# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

import numpy as np
from math import sin, cos, radians, sqrt, asin


def determine_season_ids(month):
    i_month = int(month)
    case_fnames = []
    if 3 >= i_month >= 1:
        case_fnames = ['01', '02', '03']
    elif 6 >= i_month >= 4:
        case_fnames = ['04', '05', '06']
    elif 9 >= i_month >= 7:
        case_fnames = ['07', '08', '09']
    elif 12 >= i_month >= 10:
        case_fnames = ['10', '11', '12']

    return case_fnames


def convert2julian7(year, month, day, hour, minute, second, millisecond):
    ###########################################################################

    days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days_of_month_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    days_of_month = np.asarray(days_of_month)
    days_of_month_leap = np.asarray(days_of_month_leap)

    year = int(year)
    month = int(month)
    day = int(day)
    hour = int(hour)
    minute = int(minute)
    second = int(second)

    leap = 'false'

    if (int(year/4.0) - year/4.0 == 0.):
        leap = 'true'

    if (leap == 'false'):
        julday = (np.sum(days_of_month[0:month-1]) + (day - 1)
                  + (hour*60. + minute + second/60. + millisecond
                     / (1000.*60.)) / (24.*60.))
    elif (leap == 'true'):
        julday = (np.sum(days_of_month_leap[0:month-1])+(day-1)
                  + (hour*60. + minute + second/60. + millisecond
                     / (1000.*60.)) / (24.*60.))

    return julday


###############################################################################
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1_new, lat1_new, lon2_new, lat2_new = map(
        radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2_new - lon1_new
    dlat = lat2_new - lat1_new
    a = sin(dlat/2.)**2 + cos(lat1_new) * cos(lat2_new) * sin(dlon/2.)**2
    c = 2 * asin(sqrt(a))
    # Radius of earth in kilometers is 6371
    km = 6371 * c
    return km


###############################################################################
def haversine_arr(lon_arr, lat_arr, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    npix = len(lon_arr)
    lon_arr_new = np.zeros(npix)
    lat_arr_new = np.zeros(npix)

    # convert decimal degrees to radians
    for ipix in range(0, npix):
        lon_arr_new[ipix], lat_arr_new[ipix], lon2_new, lat2_new = map(
            radians, [lon_arr[ipix], lat_arr[ipix], lon2, lat2])

    # haversine formula
    dlon = lon2_new - lon_arr_new[:]
    dlat = lat2_new - lat_arr_new[:]

    a = (np.sin(dlat[:]/2.)**2 + np.cos(lat_arr_new[:]) * np.cos(lat2_new)
         * np.sin(dlon[:]/2.)**2)
    c = 2 * np.arcsin(np.sqrt(a[:]))
    # Radius of earth in kilometers is 6371
    km = 6371 * c[:]

    return km


###############################################################################


def pres2alt(pressure):
    '''
    Determine altitude from site pressure.

    Parameters
    ----------
    pressure : numeric
        Atmospheric pressure (Pascals)

    Returns
    -------
    altitude : numeric
        Altitude in meters above sea level

    Notes
    ------
    The following assumptions are made

    ============================   ================
    Parameter                      Value
    ============================   ================
    Base pressure                  101325 Pa
    Temperature at zero altitude   288.15 K
    Gravitational acceleration     9.80665 m/s^2
    Lapse rate                     -6.5E-3 K/m
    Gas constant for air           287.053 J/(kgK)
    Relative Humidity              0%
    ============================   ================

    References
    -----------
    .. [1] "A Quick Derivation relating altitude to air pressure" from
       Portland State Aerospace Society, Version 1.03, 12/22/2004.
    '''

    alt = 44331.5 - 4946.62 * pressure ** (0.190263)

    return alt


###############################################################################
###############################################################################
###############################################################################

def alt2pres(altitude):
    '''
    Determine site pressure from altitude.

    Parameters
    ----------
    altitude : numeric
        Altitude in meters above sea level

    Returns
    -------
    pressure : numeric
        Atmospheric pressure (Pascals)

    Notes
    ------
    The following assumptions are made

    ============================   ================
    Parameter                      Value
    ============================   ================
    Base pressure                  101325 Pa
    Temperature at zero altitude   288.15 K
    Gravitational acceleration     9.80665 m/s^2
    Lapse rate                     -6.5E-3 K/m
    Gas constant for air           287.053 J/(kgK)
    Relative Humidity              0%
    ============================   ================

    References
    -----------
    .. [1] "A Quick Derivation relating altitude to air pressure" from
       Portland State Aerospace Society, Version 1.03, 12/22/2004.
    '''

    press = 100 * ((44331.514 - altitude) / 11880.516) ** (1 / 0.1902632)

    return press
