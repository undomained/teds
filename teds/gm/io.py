# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions for reading and writing geometry data."""
from astropy import units
from astropy.time import Time
from netCDF4 import Dataset
from pyquaternion import Quaternion
import datetime
import numpy as np

from .types import Geometry
from .types import Navigation
from teds.l1al1b.python.types import L1


def write_navigation(filename: str,
                     orbit_start: datetime.datetime,
                     navigation: Navigation) -> None:
    """Write navigation data to a file.

    Parameters
    ----------
    filename
        Output file path
    orbit_start
        Datetime of orbit start. Used to define the image time unit
    navigation
        Navigation data

    """
    default_fill = -32767
    nc = Dataset(filename, 'w')
    nc.title = 'Tango Carbon E2ES navigation data'
    dim_time = nc.createDimension('time', len(navigation.time))
    dim_vec = nc.createDimension('vector_elements', 3)
    dim_quat = nc.createDimension('quaternion_elements', 4)

    var = nc.createVariable('time', 'f8', dim_time, fill_value=default_fill)
    var.long_name = 'orbit vector time (seconds of day)'
    var.units = f'seconds since {orbit_start.strftime("%Y-%m-%d")}'
    var.valid_min = 0
    var.valid_max = 172800.0  # 2 x day
    var[:] = navigation.time

    var = nc.createVariable(
        'orb_pos', 'f8', (dim_time, dim_vec), fill_value=-9999999.0)
    var.long_name = 'orbit position vectors (J2000)'
    var.units = 'm'
    var.valid_min = -7200000.0
    var.valid_max = 7200000.0
    var[:] = navigation.orb_pos

    var = nc.createVariable(
        'att_quat', 'f8', (dim_time, dim_quat), fill_value=default_fill)
    var.long_name = 'Attitude quaternions (spacecraft to J2000)'
    var.units = '1'
    var.valid_min = -1.0
    var.valid_max = 1.0
    for i in range(navigation.att_quat.shape[0]):
        var[i, :] = np.roll(navigation.att_quat[i].elements, -1)

    var = nc.createVariable(
        'altitude', 'f8', (dim_time), fill_value=default_fill)
    var.long_name = 'satellite altitude'
    var.units = 'm'
    var.valid_min = 400e3
    var.valid_max = 6000e3
    var[:] = navigation.altitude

    nc.close()


def write_geometry(filename: str,
                   geometry: Geometry,
                   orbit_start: datetime.datetime,
                   l1: L1) -> None:
    """Write viewing and solar geometries and image attributes to a file.

    Parameters
    ---------
    filename
        Output NetCDF file name
    geometry
        Geometry content to be written
    orbit_start
        Datetime of orbit start. Used to define the image time unit
    l1
        L1 product containing the image timestamps

    """
    default_fill = -32767
    nc = Dataset(filename, 'w')
    nc.title = 'Tango Carbon E2ES geometry'
    n_alt, n_act = geometry.lat.shape
    dim_alt = nc.createDimension('along_track_sample', n_alt)
    dim_act = nc.createDimension('across_track_sample', n_act)
    dim_time = nc.createDimension('time', len(l1.timestamps))

    # Geometry
    var = nc.createVariable(
        'latitude', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'latitudes'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = np.rad2deg(geometry.lat)

    var = nc.createVariable(
        'longitude', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'longitudes'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = np.rad2deg(geometry.lon)

    var = nc.createVariable(
        'height', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'height from sea level'
    var.units = 'm'
    var.valid_min = -1000.0
    var.valid_max = 10000.0
    var[:] = 0.0

    var = nc.createVariable(
        'sensor_zenith', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'sensor zenith angles'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = np.rad2deg(geometry.vza)

    var = nc.createVariable(
        'sensor_azimuth', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'sensor azimuth angles'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = np.rad2deg(geometry.vaa)

    var = nc.createVariable(
        'solar_zenith', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'solar zenith angles'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = np.rad2deg(geometry.sza)

    var = nc.createVariable(
        'solar_azimuth', 'f8', (dim_alt, dim_act), fill_value=default_fill)
    var.long_name = 'solar azimuth angles'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = np.rad2deg(geometry.saa)

    # Image attributes
    var = nc.createVariable('tai_seconds', 'u4', dim_time, fill_value=0)
    var.long_name = 'detector image TAI time (seconds)'
    var.units = 'seconds since 1958-01-01 00:00:00 TAI'
    var.valid_min = np.uint(1956528000)
    var.valid_max = np.uint(2493072000)
    var[:] = l1.tai_seconds

    var = nc.createVariable('tai_subsec', 'u2', dim_time)
    var.long_name = 'detector image TAI time (subseconds)'
    var.units = '1/65536 s'
    var.valid_min = np.ushort(0)
    var.valid_max = np.ushort(65535)
    var[:] = (65535 * l1.tai_subsec).astype(np.ushort)

    var = nc.createVariable('time', 'f8', dim_time, fill_value=-32767)
    var.long_name = 'detector image time'
    var.description = 'integration start time in seconds of day'
    var.units = f'seconds since {orbit_start.strftime("%Y-%m-%d")}'
    var.valid_min = 0.0
    var.valid_max = 172800.0  # 2 x day
    var[:] = l1.timestamps

    var = nc.createVariable('day', 'f8', dim_time, fill_value=-32767)
    var.long_name = 'days since start of year'
    var.units = 's'
    var.valid_min = 0.0
    var.valid_max = 365.25
    dates = Time('1958-01-01', scale='tai') + l1.tai_seconds * units.s
    days = np.empty(len(dates))
    for i in range(len(days)):
        t = dates[i].datetime.timetuple()
        days[i] = (t.tm_yday + (t.tm_hour * 3600 + t.tm_min * 60 + t.tm_sec
                                + l1.tai_subsec[i]) / 86400)
    var[:] = days

    nc.close()


def read_navigation(filename: str) -> Navigation:
    nc = Dataset(filename)
    att_quat = np.empty(nc.dimensions['time'].size, dtype=Quaternion)
    for i in range(len(att_quat)):
        att_quat[i] = Quaternion(
            nc['att_quat'][i, 3], *nc['att_quat'][i, :3]).normalised
    return Navigation(nc['time'][:].data,
                      nc['orb_pos'][:].data,
                      att_quat,
                      nc['altitude'][:].data)


def read_geometry(filename: str) -> Geometry:
    """Read geometry (viewing and solar) from NetCDF file.

    If geolocation_data group is present in the file (e.g. in L1B
    product) then attempt to read from that. Otherwise attempt to read
    the geometry variables from the root group.

    Parameters
    ----------
    filename
        Path to a file containing geometry

    Returns
    -------
        Geometry

    """
    nc = Dataset(filename)
    if 'geolocation_data' in nc.groups:
        grp = nc['geolocation_data']
    else:
        grp = nc
    if 'height' in grp.variables:
        height = grp['height'][:].data
    else:
        height = np.zeros(grp['latitude'].shape)
    return Geometry(grp['latitude'][:].data,
                    grp['longitude'][:].data,
                    height,
                    grp['solar_zenith'][:].data,
                    grp['solar_azimuth'][:].data,
                    grp['sensor_zenith'][:].data,
                    grp['sensor_azimuth'][:].data)
