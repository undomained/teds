# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Geometry module.

Produces navigation data and geometry (viewing and solar geometries)
from orbit specification by the user and line-of-sight vectors from
the CKD.

"""
from astropy import units
from astropy.time import Time
from astropy.time import TimeDelta
from netCDF4 import Dataset
from pathlib import Path
from pyquaternion import Quaternion
from scipy.interpolate import interpn
import datetime
import numpy as np
import numpy.typing as npt

from .io import read_navigation
from .io import write_geometry
from .io import write_navigation
from .satellite import Satellite
from .sensor import Sensor
from .types import Geometry
from .types import Navigation
from teds import log
from teds.l1al1b import geolocate
from teds.l1al1b import interpolate_navigation_data
from teds.l1al1b import solar_model
from teds.l1al1b.io import read_ckd
from teds.l1al1b.types import L1
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration parameters.

    """
    if config['profile'] not in ('individual_spectra', 'orbit'):
        log.error(f'unknown geometry profile: {config["profile"]}')
        exit(1)
    ckd_path = Path(config['io_files']['ckd'])
    if not ckd_path.is_file():
        log.error(f'CKD ({ckd_path}) not found')
        exit(1)


def get_individual_spectra(config: dict, l1: L1) -> Geometry:
    """Generate GM output for individual spectra.

    First check consistencies of 'individual_spectra' input. Here we
    use a 2-dimensional (nalt, nact) data structure in an artificial
    way.

    Parameters
    ----------
    config
        Configuration dictionary
    l1
        L1 instance with timestamps. The timestamp arrays need to be
        resized for output.

    Returns
    -------
        Viewing and solar geometries.

    """
    n_alt = 1
    n_act = len(config['scene_spec']['sza'])
    log.info(f'Generating geometry for {n_act} across track locations')
    # Check_input
    for view in ('sza', 'saa', 'vza', 'vaa'):
        if n_act != len(config['scene_spec'][view]):
            log.error(
                f"input error in gm for {view} n_act ({n_act}) not equal "
                f"to {view} length ({len(config['scene_spec'][view])}).")
            exit(1)
    # Here we use the 2-dimensional data structure in an artificial way
    geometry = Geometry.from_shape((n_alt, n_act))
    # Give lon_grid and lat_grid some values such that subsequent
    # modules do not crash.
    geometry.lon[0, :] = np.deg2rad(10)
    geometry.lat[0, :] = np.deg2rad(50 + 0.0025 * np.arange(n_act))
    geometry.sza[0, :] = np.deg2rad(config['scene_spec']['sza'])
    geometry.saa[0, :] = np.deg2rad(config['scene_spec']['saa'])
    geometry.vza[0, :] = np.deg2rad(config['scene_spec']['vza'])
    geometry.vaa[0, :] = np.deg2rad(config['scene_spec']['vaa'])
    l1.navigation.time = l1.navigation.time[:n_alt]
    l1.tai_seconds = l1.tai_seconds[:n_alt]
    l1.tai_subsec = l1.tai_subsec[:n_alt]
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
    l1.navigation = Navigation.from_shape((n_time,))
    exposure_tai_start = (Time(orbit_start, scale='tai')
                          + datetime.timedelta(minutes=exposure_time_beg)
                          - Time('1958-01-01', scale='tai'))
    day_tai_start = (Time(datetime.datetime(orbit_start.year,
                                            orbit_start.month,
                                            orbit_start.day), scale='tai')
                     - Time('1958-01-01', scale='tai'))
    for i in range(n_time):
        seconds = (exposure_tai_start
                   + datetime.timedelta(seconds=i*interval)).to(units.s)
        l1.tai_seconds[i] = np.uint(seconds)
        l1.tai_subsec[i] = np.float64((
            seconds - TimeDelta(val=l1.tai_seconds[i]*units.s)).to(units.s))
        l1.navigation.time[i] = np.float64(
            (seconds - day_tai_start).to(units.s))
    return l1


def generate_attitude_quaternions(
        config: dict,
        lat_deg: npt.NDArray[np.floating],
        lon_deg: npt.NDArray[np.floating],
        vel: npt.NDArray[np.floating]) -> npt.NDArray[Quaternion]:
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
        tai_subsec = np.float64((
            tai_seconds
            - TimeDelta(val=np.uint(tai_seconds)*units.s)).to(units.s))
        # Solar model produces the J2000-ECEF quaternion so we need
        # the inverse of that.
        q_ecef_j2000 = solar_model(np.uint(tai_seconds), tai_subsec).inverse
        pos = navigation.orb_pos[i_pos]
        pos[:] = q_ecef_j2000.rotate(1e3 * pos)
        navigation.att_quat[i_pos] = q_ecef_j2000 * navigation.att_quat[i_pos]


def sensor_simulation(
        config: dict,
        sat_pos: dict,
        orbit_timestamps: npt.NDArray[np.datetime64],
        los: npt.NDArray[np.floating]) -> Geometry:
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
        config['sensor']['dwell_time']]
    return sensor.get_groundpoints(
        orbit_timestamps, dt_range[0], dt_range[1], dt_range[2], thetas)


def get_orbit(
        config: dict,
        l1: L1,
        aocs_navigation: Navigation | None) -> tuple[Navigation, Geometry]:
    """Propagate orbit and do geolocation.

    Parameters
    ----------
    config
        Configuration dictionary
    l1
        L1 product containing detector image timestamps. Required for
        interpolating navigation data from navigation to detector time
        axis.
    aocs_navigation
        Navigation data from the AOCS simulator. If not present, then
        run the SGP4 orbit simulator to generate the orbit positions
        and attitude quaternions. If present, use this instead.

    Returns
    -------
        Navigation data and viewing and solar geometries.

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
                              time_beg.day))).astype(np.floating)

    # Compute satellite position every n seconds
    if not aocs_navigation:
        log.info('Generating satellite orbit')
        satellite = Satellite(config['orbit'])
        sat_pos = satellite.compute_positions(orbit_timestamps)

    # Attitude quaternions
    if aocs_navigation:
        navigation = aocs_navigation
    else:
        log.info('Generating attitude quaternions')
        att_quat = generate_attitude_quaternions(
            config, sat_pos['lat'], sat_pos['lon'], sat_pos['v'])
        # Navigation data with orbit positions in ECEF. These will be
        # converted to J2000 in convert_to_j2000.
        navigation = Navigation(time=orb_time_day,
                                orb_pos=sat_pos['p'],   # ECEF here
                                att_quat=att_quat,   # SC-to-ECEF here
                                altitude=1e3 * sat_pos['height'])  # m

    # Do geolocation using the orbit positions and line-of-sight (LOS)
    # vectors from the CKD. Default is to derive the orbit positions
    # and quaternions in the J2000 frame and perform
    # geolocation. Otherwise use the fallback Python implementation.
    ckd = read_ckd(config['io_files']['ckd'])
    log.info('Geolocation')
    if not config['use_python_geolocation']:
        if not aocs_navigation:
            # Convert orbit positions and quaternions from ECEF to
            # J2000. AOCS data is already in J2000.
            convert_to_j2000(orbit_timestamps, navigation)
        l1.navigation = interpolate_navigation_data(
            navigation, l1.navigation.time)
        geometry = geolocate(
            l1, ckd.swath.line_of_sights, config['io_files']['dem'])
    else:
        # Configure sensors and compute ground pixel information
        geometry = sensor_simulation(
            config, sat_pos, orbit_timestamps, ckd.swath.line_of_sights)
    return navigation, geometry


def extend_geometry(geometry: Geometry,
                    margin_alt: int,
                    margin_act: int,
                    density_alt: int,
                    density_act: int) -> Geometry:
    """Extend geometry outside the target box.

    Parameters
    ----------
    geometry
        Input geometry. This is the same as in the L1B product.
    margin_alt
        Number of grid points, in units of the input grid, outside the
        target box in the ALT direction. This applies separately to
        both ends of the geometry, i.e. the total number of out-of-box
        points is 2 * margin_alt.
    margin_act
        Similarly, number of added grid points in the ACT direction.
    density_alt
        Number of extended grid points as a multiple of the input grid
        in the ALT direction.
    density_act
        Number of extended grid points as a multiple of the input grid
        in the ACT direction.

    Returns
    -------
        Extended and densified geometry

    """
    n_rows, n_cols = geometry.lat.shape
    alt_in = np.arange(0, n_rows, 1.0)
    act_in = np.arange(0, n_cols, 1.0)
    grid_out = np.meshgrid(
        np.arange(-margin_alt, n_rows + margin_alt, 1.0 / density_alt),
        np.arange(-margin_act, n_cols + margin_act, 1.0 / density_act),
        indexing='ij')

    def interp(data: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Interpolate to a target (dense) grid."""
        return interpn((alt_in, act_in),
                       data,
                       grid_out,
                       method='linear',
                       bounds_error=False,
                       fill_value=None)

    def interp_angles(data: npt.NDArray[np.floating]) -> (
            npt.NDArray[np.floating]):
        """For angles, need to interpolate sin and cos separately."""
        return np.arctan2(interp(np.sin(data)), interp(np.cos(data)))

    return Geometry(interp(geometry.lat),
                    interp_angles(geometry.lon),
                    interp(geometry.height),
                    interp(geometry.sza),
                    interp_angles(geometry.saa),
                    interp(geometry.vza),
                    interp_angles(geometry.vaa))


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
    print(flush=True)

    config = merge_config_with_default(config_user, 'teds.gm')
    check_config(config)

    # If present, use AOCS navigation data. Otherwise ignore it and
    # simulate orbit using the SPG4 routines instead.
    aocs_navigation = None
    if config['io_files']['aocs_navigation']:
        log.info('Reading AOCS generated quaternions')
        aocs_navigation = read_navigation(
            config['io_files']['aocs_navigation'])
        config['orbit']['epoch'] = datetime.datetime.fromisoformat(Dataset(
            config['io_files']['aocs_navigation']).epoch)

    # Generate image timestamps which are essentially part of geometry
    # information. Without timestamps it is not clear what the
    # retrieved latitudes and longitudes refer to.
    log.info('Generating detector image timestamps')
    l1 = gen_image_timestamps(config['orbit']['epoch'],
                              config['sensor']['start_time'],
                              config['sensor']['end_time'],
                              config['sensor']['dwell_time'])

    if config['profile'] == 'individual_spectra':
        geometry = get_individual_spectra(config, l1)
        # Overwrite these settings. Extended geometry does not apply here.
        config['extend_geometry']['margin_alt'] = 0
        config['extend_geometry']['margin_act'] = 0
        config['extend_geometry']['density_alt'] = 1
        config['extend_geometry']['density_act'] = 1
    elif config['profile'] == 'orbit':
        navigation, geometry = get_orbit(config, l1, aocs_navigation)
    else:
        log.error(f'unknown profile: {config["profile"]}')
        exit(1)

    # Extended geometry may be used in the SGM
    geometry_ext = extend_geometry(geometry,
                                   config['extend_geometry']['margin_alt'],
                                   config['extend_geometry']['margin_act'],
                                   config['extend_geometry']['density_alt'],
                                   config['extend_geometry']['density_act'])

    # Write output data
    log.info('Writing output')
    write_geometry(config['io_files']['geometry'],
                   geometry,
                   geometry_ext,
                   config['orbit']['epoch'],
                   l1.navigation.time,
                   l1.tai_seconds,
                   l1.tai_subsec,)
    if config['profile'] == 'orbit' and not config['use_python_geolocation']:
        write_navigation(config['io_files']['navigation'],
                         config['orbit']['epoch'],
                         navigation)

    # If this is shown then the simulation ran successfully
    print_heading('Success')
