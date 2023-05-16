# library GM functions

import math
import sys
import numpy as np
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


