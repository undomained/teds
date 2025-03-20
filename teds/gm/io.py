# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions for reading and writing geometry data."""
from astropy import units
from astropy.time import Time
from netCDF4 import Dataset
from netCDF4 import Group
from pyquaternion import Quaternion
import datetime
import numpy as np
import numpy.typing as npt

from .types import Geometry
from .types import Navigation

default_fill = -32767


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


def nc_write_lat(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('latitude',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'latitude at bin locations'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = data


def nc_write_lon(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('longitude',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'longitudes at bin locations'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = data


def nc_write_height(nc: Dataset | Group,
                    data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('height',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'height from sea level at bin locations'
    var.units = 'm'
    var.valid_min = -1000.0
    var.valid_max = 10000.0
    var[:] = data


def nc_write_vza(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('sensor_zenith',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'sensor zenith angle at bin locations'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = data


def nc_write_vaa(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('sensor_azimuth',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'sensor azimuth angle at bin locations'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = data


def nc_write_sza(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('solar_zenith',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'solar zenith angle at bin locations'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = data


def nc_write_saa(nc: Dataset | Group, data: npt.NDArray[np.float64]) -> None:
    var = nc.createVariable('solar_azimuth',
                            'f8',
                            ('along_track_sample', 'across_track_sample'),
                            fill_value=default_fill)
    var.long_name = 'solar azimuth angle at bin locations'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = data


def nc_write_geometry(nc: Dataset | Group,
                      geometry: Geometry,
                      convert_rad2deg: bool = True) -> None:
    """Write geometry to a NetCDF group.

    This wraps the other nc_write_* functions and can be used if the
    geometry is complete (all variables have content).

    Parameters
    ----------
    nc
        NetCDF group (Dataset if root group)
    geometry
        Object representing latitudes, longitudes, and solar and
        viewing geometries. All angles are in radians.
    convert_rad2deg
        Whether to convert to degrees

    """
    if convert_rad2deg:
        nc_write_lat(nc, np.rad2deg(geometry.lat))
        nc_write_lon(nc, np.rad2deg(geometry.lon))
        nc_write_vza(nc, np.rad2deg(geometry.vza))
        nc_write_vaa(nc, np.rad2deg(geometry.vaa))
        nc_write_sza(nc, np.rad2deg(geometry.sza))
        nc_write_saa(nc, np.rad2deg(geometry.saa))
    else:
        nc_write_lat(nc, geometry.lat)
        nc_write_lon(nc, geometry.lon)
        nc_write_vza(nc, geometry.vza)
        nc_write_vaa(nc, geometry.vaa)
        nc_write_sza(nc, geometry.sza)
        nc_write_saa(nc, geometry.saa)
    nc_write_height(nc, np.zeros(geometry.lat.shape))


def write_geometry(filename: str,
                   geometry: Geometry,
                   geometry_ext: Geometry,
                   orbit_start: datetime.datetime,
                   timestamps: npt.NDArray[np.float64],
                   tai_seconds: npt.NDArray[np.uint],
                   tai_subsec: npt.NDArray[np.float64]) -> None:
    """Write viewing and solar geometries and image attributes to a file.

    Parameters
    ---------
    filename
        Output NetCDF file name
    geometry
        Geometry content to be written
    orbit_start
        Datetime of orbit start. Used to define the image time unit
    timestamps
        Detector image timestamps
    tai_seconds
        Number of TAI seconds since the beginning of epoch
        corresponding to detector timestmps (integer part).
    tai_subsec
        Fractional part of TAI seconds since the beginning of epoch.

    """
    nc = Dataset(filename, 'w')
    nc.title = 'Tango Carbon E2ES geometry'
    n_alt, n_act = geometry.lat.shape
    dim_alt = nc.createDimension('along_track_sample', n_alt)
    nc.createDimension('across_track_sample', n_act)

    # Normal geometry
    nc_write_geometry(nc, geometry)

    # Extended geometry (for SGM)
    grp = nc.createGroup('extended')
    grp.createDimension('along_track_sample', geometry_ext.lat.shape[0])
    grp.createDimension('across_track_sample', geometry_ext.lat.shape[1])
    nc_write_geometry(grp, geometry_ext)

    # Image attributes
    var = nc.createVariable('tai_seconds', 'u4', dim_alt, fill_value=0)
    var.long_name = 'detector image TAI time (seconds)'
    var.units = 'seconds since 1958-01-01 00:00:00 TAI'
    var.valid_min = np.uint(1956528000)
    var.valid_max = np.uint(2493072000)
    var[:] = tai_seconds

    var = nc.createVariable('tai_subsec', 'u2', dim_alt)
    var.long_name = 'detector image TAI time (subseconds)'
    var.units = '1/65536 s'
    var.valid_min = np.ushort(0)
    var.valid_max = np.ushort(65535)
    var[:] = (65535 * tai_subsec).astype(np.ushort)

    var = nc.createVariable('time', 'f8', dim_alt, fill_value=-32767)
    var.long_name = 'detector image time'
    var.description = 'integration start time in seconds of day'
    var.units = f'seconds since {orbit_start.strftime("%Y-%m-%d")}'
    var.valid_min = 0.0
    var.valid_max = 172800.0  # 2 x day
    var[:] = timestamps

    var = nc.createVariable('day', 'f8', dim_alt, fill_value=-32767)
    var.long_name = 'days since start of year'
    var.units = 's'
    var.valid_min = 0.0
    var.valid_max = 365.25
    dates = Time('1958-01-01', scale='tai') + tai_seconds * units.s
    days = np.empty(len(dates))
    for i in range(len(days)):
        t = dates[i].datetime.timetuple()
        days[i] = (t.tm_yday + (t.tm_hour * 3600 + t.tm_min * 60 + t.tm_sec
                                + tai_subsec[i]) / 86400)
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


def read_geometry(filename: str, group: str | None = None) -> Geometry:
    """Read geometry (viewing and solar) from NetCDF file.

    If geolocation_data group is present in the file (e.g. in L1B
    product) then attempt to read from that. Otherwise attempt to read
    the geometry variables from the root group.

    Parameters
    ----------
    filename
        Path to a file containing geometry
    group
        If present, read from this NetCDF group. Default is the root
        group.

    Returns
    -------
        Geometry

    """
    nc = Dataset(filename)
    if group:
        grp = nc[group]
    elif 'geolocation_data' in nc.groups:
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
