#==============================================================================
#   Tools for geometry module
#   This source code is licensed under the 3-clause BSD license found in
#   the LICENSE file in the root directory of this project.
#==============================================================================

import math
import sys
import numpy as np
import numpy.typing as npt
from pyproj import Geod
from scipy.io import netcdf
import netCDF4 as nc


def sunpos(year, month, day, hour, minute, second, timezone, latitude, longitude, refraction):

    # @author: John Clark Craig
    # https://levelup.gitconnected.com/python-sun-position-for-solar-energy-and-research-7a4ead801777

    # Math typing shortcuts
    rad, deg = math.radians, math.degrees
    sin, cos, tan = math.sin, math.cos, math.tan
    asin, atan2 = math.asin, math.atan2

# Decimal hour of the day at Greenwich
    greenwichtime = hour - timezone + minute / 60 + second / 3600

# Days from J2000, accurate from 1901 to 2099
    daynum = (
        367 * year
        - 7 * (year + (month + 9) // 12) // 4
        + 275 * month // 9
        + day
        - 730531.5
        + greenwichtime / 24
    )
# Mean longitude of the sun
    mean_long = daynum * 0.01720279239 + 4.894967873
# Mean anomaly of the Sun
    mean_anom = daynum * 0.01720197034 + 6.240040768
# Ecliptic longitude of the sun
    eclip_long = (
        mean_long
        + 0.03342305518 * sin(mean_anom)
        + 0.0003490658504 * sin(2 * mean_anom)
    )
# Obliquity of the ecliptic
    obliquity = 0.4090877234 - 0.000000006981317008 * daynum
# Right ascension of the sun
    rasc = atan2(cos(obliquity) * sin(eclip_long), cos(eclip_long))
# Declination of the sun
    decl = asin(sin(obliquity) * sin(eclip_long))

# Local elevation of the sun
    sza = []
    saz = []
    for lat, lon in zip(latitude, longitude):
        rlat = rad(lat)
        rlon = rad(lon)

# Local sidereal time
        sidereal = 4.894961213 + 6.300388099 * daynum + rlon

# Hour angle of the sun
        hour_ang = sidereal - rasc

        elevation = asin(sin(decl) * sin(rlat) + cos(decl) * cos(rlat) * cos(hour_ang))

# Local azimuth of the sun
        azimuth = atan2(
            -cos(decl) * cos(rlat) * sin(hour_ang),
            sin(decl) - sin(rlat) * sin(elevation),
        )
# Convert azimuth and elevation to degrees
        azimuth = into_range(deg(azimuth), 0, 360)
        elevation = into_range(deg(elevation), -180, 180)
# Refraction correction (optional)
        if refraction:
            targ = rad((elevation + (10.3 / (elevation + 5.11))))
            elevation += (1.02 / tan(targ)) / 60
# Return azimuth and elevation in degrees
        sza = np.append(sza, 90.-elevation)
        saz = np.append(saz, azimuth)
    return (saz, sza)


def into_range(x, range_min, range_max):
    shiftedx = x - range_min
    delta = range_max - range_min
    return (((shiftedx % delta) + delta) % delta) + range_min


def trans_lat_lon(lat_ref, lon_ref, when, spat_grid, nalt, nact):

    year, month, day, hour, minute, timezone = when
    salt, sact = spat_grid

    lat_grid = np.empty([nalt, nact])
    lon_grid = np.empty([nalt, nact])
    saz = np.empty([nalt, nact])
    sza = np.empty([nalt, nact])

    # sact is the sampling distance from subsatellite point to the ground pixel in across track direction.
    # salt is the along track sampling at the ACT pixel. Due to sphericity, this varies over the swath.
    # salt and sact are both of dimension nact.

    # latitude and longitude reference coordinates, needs to be dublicated for efficient g.fwd simulation
    lon_ref = np.zeros(nact)+lon_ref
    lat_ref = np.zeros(nact)+lat_ref

    # some settings for the sza and saz calculation
    second = 0  # do not discreminate seconds....
    refraction = True

    # g = Geod(ellps='clrk66') # Use Clarke 1866 ellipsoid.
    g = Geod(ellps='WGS84')  # The World Geodetic System 1984 (WGS84) reference

    # we assume that the reference point is located at scan line nalt/2
    # thus, first scan line is at a along track distance -nalt/2 * salt
    # so each of a consecuative scan line is at -nalt/2 * salt + ialt*salt = ialt-nalt/2

    index = -np.int0(nalt/2 + 0.5)

    for ialt in range(0, nalt):
        alt_distance = index*salt
        distance = np.sqrt(sact**2 + (alt_distance)**2)
        azi_fwd = np.arctan2(sact, alt_distance)/np.pi*180.
        lon_grid[ialt, :], lat_grid[ialt, :], azi_fwd = g.fwd(lon_ref, lat_ref, azi_fwd, distance)
        saz[ialt, :], sza[ialt, :] = sunpos(year, month, day, hour, minute, second,
                                            timezone, lat_grid[ialt], lon_grid[ialt], refraction)
        index += 1

    return(sza, saz, lat_grid, lon_grid)


def haversine_formula(lat1, lon1, lat2, lon2, rad_earth):
  #  lat1 and lat2: degree
  #  lon1 and lon2: degree
  #  rad_earth: Earth radius in m

    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = np.sin(dlat/2.)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon/2.)**2
    # Haversine formula
    c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))
    return rad_earth * c

def vincenty(lat1: npt.NDArray[np.float64],
             lat2: npt.NDArray[np.float64],
             lon1: npt.NDArray[np.float64],
             lon2: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    a = 6378137.0
    b = 6356752.314245179499
    f = 3.3528106647474804385e-3
    U1 = np.arctan((1 - f) * np.tan(lat1))
    U2 = np.arctan((1 - f) * np.tan(lat2))
    sU1 = np.sin(U1)
    sU2 = np.sin(U2)
    cU1 = np.cos(U1)
    cU2 = np.cos(U2)
    L = lon2 - lon1
    cos_alpha2 = 0
    sigma = 0
    sin_sigma = 0
    cos_sigma = 0
    cos_2sig = 0
    _lambda = L
    _lambda_prev = _lambda
    max_iter = 10
    for i_iter in range(max_iter):
        slam = np.sin(_lambda)
        clam = np.cos(_lambda)
        cos_U2_t1 = cU2 * slam
        cos_U2_t2 = cU1 * sU2 - sU1 * cU2 * clam
        sin_sigma = np.sqrt(cos_U2_t1 * cos_U2_t1 + cos_U2_t2 * cos_U2_t2)
        cos_sigma = sU1 * sU2 + cU1 * cU2 * clam
        sigma = np.arctan2(sin_sigma, cos_sigma)
        sin_alpha = np.where(np.abs(sin_sigma < 1e-30),
                             cU1,
                             cU1 * cU2 * slam / sin_sigma)
 #       sin_alpha = cU1 * cU2 * slam / sin_sigma
        cos_alpha2 = 1 - sin_alpha * sin_alpha
        cos_2sig = cos_sigma - 2 * sU1 * sU2 / cos_alpha2
        C = f / 16 * cos_alpha2 * (4.0 + f * (4 - 3 * cos_alpha2))
        _lambda = (L + (1 - C) * f * sin_alpha
                   * (sigma + C * sin_sigma
                      * (cos_2sig + C * cos_sigma * (-1 + 2 * cos_2sig))))
        if (abs(_lambda - _lambda_prev).max() < 1e-12):
            break
        _lambda_prev = _lambda
    u2 = cos_alpha2 * (a * a - b * b) / (b * b)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    cos_2sig2 = cos_2sig * cos_2sig
    D_sigma = (B * sin_sigma
               * (cos_2sig
                  + 1.0 / 4 * B
                  * (cos_sigma * (-1 + 2 * cos_2sig2)
                     - B / 6 * cos_2sig
                     * (-3 + 4 * sin_sigma * sin_sigma)
                     * (-3 + 4 * cos_2sig2))))
    return b * A * (sigma - D_sigma)

