# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Geometry module.

Produces navigation data and geometry (viewing and solar geometries)
from orbit specification by the user and line-of-sight vectors from
the CKD.

"""
from astropy import units
from astropy.time import Time
from pathlib import Path
from pyquaternion import Quaternion
from scipy.interpolate import CubicSpline
import datetime
import numpy as np
import numpy.typing as npt

from .io import write_geometry
from .io import write_image_attributes
from .io import write_navigation
from .satellite import Satellite
from .sensor import Sensor
from .types import Navigation
from teds import log
from teds.l1al1b import geolocate
from teds.l1al1b import solar_model
from teds.l1al1b.python.io import read_ckd
from teds.l1al1b.python.types import Geometry
from teds.l1al1b.python.types import L1
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config_file
        Path of YAML configuration file.

    """
    if config['profile'] not in ('individual_spectra', 'orbit'):
        log.error(f'unknown geometry profile: {config["profile"]}')
        exit(1)
    ckd_path = Path(config['io_files']['ckd'])
    if not ckd_path.is_file():
        log.error(f'CKD ({ckd_path}) not found')
        exit(1)


def get_individual_spectra(config: dict) -> Geometry:
    """Generate GM output for individual spectra.

    First check consistencies of 'individual_spectra' input. Here we
    use a 2-dimensional (nalt, nact) data structure in an artificial
    way.

    Parameters
    ----------
    config
        Configuration dictionary

    Returns
    -------
        Viewing and solar geometries.

    """
    nalt = 1
    nact = len(config['scene_spec']['sza'])
    log.info(f'Generating geometry for {nact} across track locations')
    # Check_input
    for view in ('sza', 'saa', 'vza', 'vaa'):
        if nact != len(config['scene_spec'][view]):
            log.error(f"input error in gm for {view} nact ({nact}) not equal "
                      f"to {view} length ({len(config['scene_spec'][view])}).")
            exit(1)
    # Here we use the 2-dimensional data structure in an artificial way
    geometry = Geometry.from_shape((nalt, nact))
    # Give lon_grid and lat_grid some values such that subsequent
    # modules do not crash.
    geometry.lon[0, :] = np.deg2rad(10)
    geometry.lat[0, :] = np.deg2rad(50 + 0.0025 * np.arange(nact))
    geometry.sza[0, :] = np.deg2rad(config['scene_spec']['sza'])
    geometry.saa[0, :] = np.deg2rad(config['scene_spec']['saa'])
    geometry.vza[0, :] = np.deg2rad(config['scene_spec']['vza'])
    geometry.vaa[0, :] = np.deg2rad(config['scene_spec']['vaa'])
    return geometry


def gen_orbit_timestamps(dt_begin: datetime.datetime,
                         duration: float) -> npt.NDArray[np.datetime64]:
    """Generate orbit timestamps.

    Parameters
    ----------
    datetime_begin
        Datetime of the first orbit position
    duration
        orbit duration in hours

    Returns
    -------
        List of orbit timestamps.

    """
    dt_end = dt_begin + datetime.timedelta(hours=duration)
    dt_interval = 10.0
    return np.arange(dt_begin, dt_end, datetime.timedelta(seconds=dt_interval))


def gen_image_timestamps(orbit_start: datetime.datetime,
                         exposure_time_beg: float,
                         exposure_time_end: float,
                         interval: float) -> L1:
    """Generate image timestamps.

    Parameters
    ----------
    orbit start
        Datetime of the beginning of orbit.
    exposure_time_beg
        Time of first detector exposure since the beginning of orbit,
        in minutes.
    exposure_time_end
        Time of last detector exposure since the beginning of orbit,
        in minutes.
    interval
        Exposure time in seconds.

    Returns
    -------
        L1 product which is mostly empty but contains number of TAI
        seconds since the beginning of TAI epoch (split into whole and
        fractional parts) and image time in units of seconds since the
        beginning of day.

    """
    n_time = int(60 * (exposure_time_end - exposure_time_beg) / interval)
    # L1 products, only used for storing the detector image times
    l1 = L1.from_empty()
    l1.tai_seconds = np.empty(n_time, dtype=np.uint32)
    l1.tai_subsec = np.empty(n_time)
    l1.timestamps = np.empty(n_time, dtype=np.float64)
    exposure_tai_start = (Time(orbit_start, scale='tai')
                          + datetime.timedelta(minutes=exposure_time_beg)
                          - Time('1958-01-01', scale='tai'))
    day_tai_start = (Time(datetime.datetime(orbit_start.year,
                                            orbit_start.month,
                                            orbit_start.day), scale='tai')
                     - Time('1958-01-01', scale='tai'))
    # Use hours since 1958 to offset the TAI time so that the
    # subsecond component can be calculated with better numerical
    # precision.
    hours_since_epoch = int(np.float64(day_tai_start.to(units.h)))
    for i in range(n_time):
        tai_seconds = (exposure_tai_start
                       + datetime.timedelta(seconds=i*interval)).to(units.s)
        l1.tai_seconds[i] = np.uint(tai_seconds)
        # Get the fractional part of TAI seconds. In order to improve
        # numerical precision, reduce the total number of seconds by
        # offsetting with hours since epoch.
        tai_seconds_reduced = (
            exposure_tai_start
            + datetime.timedelta(seconds=i*interval)
            - datetime.timedelta(hours=hours_since_epoch)).to(units.s)
        l1.tai_subsec[i] = np.float64(tai_seconds_reduced) % 1
        l1.timestamps[i] = np.float64((tai_seconds
                                       - day_tai_start).to(units.s))
    return l1


def generate_attitude_quaternions(
        config: dict,
        lat_deg: npt.NDArray[np.float64],
        lon_deg: npt.NDArray[np.float64],
        vel: npt.NDArray[np.float64]) -> npt.NDArray[Quaternion]:
    """Generate nominal attitude quaternions.

    The quaternions correspond to a rotation of the nadir
    line-of-sight (LOS) vector from the spacecraft (SC) reference
    frame to the Earth-centered Earth-fixed (ECEF) frame. The strategy
    is to first find the ECEF-to-SC quaternion, which is the other way
    around because it is easier, and then conjugate that.

    Find first the rotation that aligns the z-axes and then the
    rotation that aligns the x-axes. Rotation to align one vector with
    another does not have a unique solution but the shortest rotation
    is given by
      qz = (|v1||v2| + v1 v2) + v1 x v2,
    where the first term is the scalar part of the quaternion. With v1
    the nadir vector in ECEF and v2 the nadir vector in SC frame, the
    rotation to align z-axes is
      qz = 1 + Nzz + Nz x (0, 0, 1).
    Apply this rotation to the x-axis in ECEF,
      Nx^(z) = qz Nx qz^-1,
    and find the rotation to align the x-axes:
      qx = 1 + Nxx^(z) + Nx^(z) x (1, 0, 0).
    The required SC-to-ECEF rotation is the inverse of the two
    rotations:
      q^SC-to-ECEF = (qx qz)^-1.

    Parameters
    ----------
    config
        Configuration dictionary
    lat_deg
        Orbit ground track latitudes in degrees
    lon_deg
        Orbit ground track longitudes in degrees
    vel
        Spacecraft velocities

    Returns
    -------
        Attitude quaternions in the form of a Numpy array

    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    # The original satellite velocity vectors. Might not be exactly
    # orthogonal to the nadir direction.
    sat_x_orig = vel / np.linalg.norm(vel, axis=1)[:, None]
    # The nadir vector. Minus sign because we are looking from the
    # spacecraft down to Earth.
    sat_z = np.empty((len(lat), 3))
    sat_z[:, 0] = -np.cos(lat) * np.cos(lon)
    sat_z[:, 1] = -np.cos(lat) * np.sin(lon)
    sat_z[:, 2] = -np.sin(lat)
    # Ensure the velocity and nadir directions are orthogonal
    cross = np.cross(sat_x_orig, sat_z)
    sat_x = np.cross(sat_z, cross)
    att_quat = np.empty(len(lat), dtype=Quaternion)
    for i_att in range(len(att_quat)):
        u_z = [0, 0, 1]
        vec = np.cross(sat_z[i_att, :], u_z)
        # Scalar part is sqrt(u_z^2 * sat_z^2) + u_z * sat_z
        scalar = 1 + sat_z[i_att, 2]
        q_z = Quaternion(scalar, *vec).normalised
        sat_x_rotated_z = q_z.rotate(sat_x[i_att, :])
        u_x = [1, 0, 0]
        vec = np.cross(sat_x_rotated_z, u_x)
        scalar = 1 + np.dot(u_x, sat_x_rotated_z)
        q_x = Quaternion(scalar, *vec).normalised
        q = q_x * q_z
        # Take the inverse because we need to rotate in the other
        # direction.
        att_quat[i_att] = q.inverse.elements
    return att_quat


def convert_to_j2000(orbit_timestamps: npt.NDArray[np.datetime64],
                     navigation: Navigation) -> None:
    """Convert orbit positions and quaternions from ECEF to J2000.

    Parameters
    ----------
    l1b_lib
        L1B processor C++ Library
    orbit_timestamps
        Orbit timestamps
    Navigation
        Orbit positions and nominal attitude quaternions in ECEF

    """
    for i_pos in range(len(navigation.orb_pos)):
        tai_seconds = (Time(orbit_timestamps[i_pos], scale='tai')
                       - Time('1958-01-01', scale='tai')).to(units.s)
        tai_subsec = np.float64(np.float128(tai_seconds) % 1)
        # Solar model produces the J2000-ECEF quaternion so we need
        # the inverse of that.
        q_ecef_j2000 = solar_model(np.uint(tai_seconds), tai_subsec).inverse
        pos = navigation.orb_pos[i_pos]
        pos[:] = q_ecef_j2000.rotate(1e3 * pos)
        navigation.att_quat[i_pos] = q_ecef_j2000 * navigation.att_quat[i_pos]


def interpolate_navigation_data(navigation: Navigation, l1: L1) -> None:
    """Interpolate navigation data from orbit times to detector image
    times.

    """
    n_alt = len(l1.timestamps)
    l1.orb_pos = np.empty((n_alt, 3))
    l1.att_quat = np.empty(n_alt, dtype=Quaternion)
    # Interpolate orbit positions
    for i_dir in range(3):
        s = CubicSpline(navigation.time, navigation.orb_pos[:, i_dir])
        l1.orb_pos[:, i_dir] = s(l1.timestamps)
    # Interpolate quaternions
    indices = np.searchsorted(
        navigation.time.astype(np.float64), l1.timestamps)
    for i_alt in range(n_alt):
        idx_lo = indices[i_alt] - 1
        idx_hi = indices[i_alt]
        idx_delta = ((l1.timestamps[i_alt] - navigation.time[idx_lo])
                     / (navigation.time[idx_hi] - navigation.time[idx_lo]))
        q0 = navigation.att_quat[idx_lo]
        q1 = navigation.att_quat[idx_hi]
        l1.att_quat[i_alt] = Quaternion.slerp(q0, q1, amount=idx_delta)


def sensor_simulation(
        config: dict,
        sat_pos: dict,
        orbit_timestamps: npt.NDArray[np.datetime64],
        los: npt.NDArray[np.float64]) -> Geometry:
    """Propogate sensor."""
    thetas = np.rad2deg(np.arctan(los[:, 1] / los[:, 2]))
    # Make and propage the sensor
    sensor = Sensor(thetas)
    sensor.compute_groundpoints(sat_pos,
                                pitch=config['sensor']['pitch'],
                                yaw=config['sensor']['yaw'],
                                roll=config['sensor']['roll'])
    log.info('Computing ground pixels:')
    # Time range of observations and interval
    dt_range = [
        orbit_timestamps[0].astype(datetime.datetime)
        + datetime.timedelta(minutes=config['sensor']['start_time']),
        orbit_timestamps[0].astype(datetime.datetime)
        + datetime.timedelta(minutes=config['sensor']['end_time']),
        config['sensor']['integration_time']]
    return sensor.get_groundpoints(
        orbit_timestamps, dt_range[0], dt_range[1], dt_range[2], thetas)


def get_orbit(config: dict) -> tuple[Navigation, Geometry, L1]:
    """Propagate orbit and do geolocation.

    Parameters
    ----------
    config
        Configuration dictionary

    Returns
    -------
        Navigation data, viewing and solar geometries, and detector
        image timestamps. The latter are essentially part of geometry
        information. Without the timestamps it is not clear what the
        retrieved latitudes and longitudes refer to.

    """
    # Orbit times as datetime objects
    log.info('Generating orbit timestamps')
    orbit_timestamps = gen_orbit_timestamps(
        config['orbit']['epoch'], config['orbit']['propagation_duration'])
    # Orbit times in seconds since beginning of day
    time_beg = orbit_timestamps[0].astype(datetime.datetime)
    orb_time_day = 1e-6 * (
        orbit_timestamps - np.datetime64(
            datetime.datetime(time_beg.year,
                              time_beg.month,
                              time_beg.day))).astype(np.float64)

    # Compute satellite position every n seconds
    log.info('Generating satellite orbit')
    satellite = Satellite(config['orbit'])
    sat_pos = satellite.compute_positions(orbit_timestamps)

    # Generate image timestamps
    log.info('Generating detector image timestamps')
    l1 = gen_image_timestamps(config['orbit']['epoch'],
                              config['sensor']['start_time'],
                              config['sensor']['end_time'],
                              config['sensor']['integration_time'])

    # Attitude quaternions
    log.info('Generating attitude quaternions')
    att_quat = generate_attitude_quaternions(
        config, sat_pos['lat'], sat_pos['lon'], sat_pos['v'])

    # Navigation data with orbit positions in ECEF. These will be
    # converted to J2000 in convert_to_j2000.
    navigation = Navigation(time=orb_time_day,
                            orb_pos=sat_pos['p'],   # ECEF for now
                            att_quat=att_quat,   # SC-to-ECEF for now
                            altitude=1e3 * sat_pos['height'])  # m

    # Do geolocation using the orbit positions and line-of-sight (LOS)
    # vectors from the CKD. Default is to derive the orbit positions
    # and quaternions in the J2000 frame and perform
    # geolocation. Otherwise use the fallback Python implementation.
    ckd = read_ckd(config['io_files']['ckd'])
    log.info('Geolocation')
    if not config['use_python_geolocation']:
        # Convert orbit positions and quaternions from ECEF to J2000
        convert_to_j2000(orbit_timestamps, navigation)
        interpolate_navigation_data(navigation, l1)
        geometry = geolocate(
            l1, ckd.swath.line_of_sights, config['io_files']['dem'])
    else:
        # Configure sensors and compute ground pixel information
        geometry = sensor_simulation(
            config, sat_pos, orbit_timestamps, ckd.swath.line_of_sights)
    return navigation, geometry, l1


def geometry_module(config_user: dict | None = None) -> None:
    """Generate viewing and solar geometries from orbit specification.

    If run without an argument, print the list of all settings and
    exit.

    Parameters
    ----------
    config_user
        Configuration dictionary directly from file, as given by the
        user, to be expanded with default values for parameters not
        specified by the user.

    """
    print_heading('Tango geometry module', empty_line=False)
    print_system_info()
    print()

    config = merge_config_with_default(config_user, 'teds.gm')
    check_config(config)

    if config['profile'] == 'individual_spectra':
        geometry = get_individual_spectra(config)
    elif config['profile'] == 'orbit':
        navigation, geometry, l1 = get_orbit(config)
    else:
        log.error(f'unknown profile: {config["profile"]}')
        exit(1)

    # Write output data
    write_geometry(config['io_files']['geometry'], geometry)
    if config['profile'] == 'orbit' and not config['use_python_geolocation']:
        write_image_attributes(
            config['io_files']['geometry'], config['orbit']['epoch'], l1)
        write_navigation(config['io_files']['navigation'],
                         config['orbit']['epoch'],
                         navigation)

    # If this is shown then the simulation ran successfully
    print_heading('Success')
