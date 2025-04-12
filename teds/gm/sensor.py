# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Class for simulating a Tango instrument (Carbon/Nitro) orientation.

This is the fallback implementation if the C++ geolocation functions
are not available.

"""
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_sun
from astropy.time import Time
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
from sgp4.earth_gravity import wgs84
import datetime
import numpy as np
import numpy.typing as npt

from .types import Geometry
from teds import log


class Sensor:
    def __init__(self, thetas: npt.NDArray[np.floating]) -> None:
        self.thetas = thetas
        self.n_thetas = self.thetas.shape[0]

    @staticmethod
    def trans_yaw_pitch_roll(y: float,
                             p: float,
                             r: float,
                             degrees: bool = True) -> npt.NDArray[np.floating]:
        """Transformation matrix.

        Ordering is: yaw, pitch, roll.

        """
        if degrees:
            y = np.deg2rad(y)
            p = np.deg2rad(p)
            r = np.deg2rad(r)
        sin_y = np.sin(y)
        cos_y = np.cos(y)
        sin_p = np.sin(p)
        cos_p = np.cos(p)
        sin_r = np.sin(r)
        cos_r = np.cos(r)
        rot_mat = np.array([
            [cos_y*cos_p,
             cos_y*sin_p*sin_r - sin_y*cos_r,
             cos_y*sin_p*cos_r + sin_y*sin_r],
            [sin_y*cos_p,
             sin_y*sin_p*sin_r + cos_y*cos_r,
             sin_y*sin_p*cos_r - cos_y*sin_r],
            [-sin_p, cos_p*sin_r, cos_p*cos_r]])
        return rot_mat

    @staticmethod
    def scale_axis(x: npt.NDArray[np.floating]) -> tuple[
            npt.NDArray[np.floating], np.floating, np.floating]:
        """Scale axis from -1 to 1."""
        offset = np.mean(x)
        scale_factor = 0.5 * (np.max(x) - np.min(x))
        return (x - offset) / scale_factor, offset, scale_factor

    @staticmethod
    def apply_scale_axis(
            x: npt.NDArray[np.floating],
            scaled_axis: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Apply scale factor and offset."""
        return (x - scaled_axis[1]) / scaled_axis[2]

    @staticmethod
    def line_sphere_intersect(
            p_sat: npt.NDArray[np.floating],
            p: npt.NDArray[np.floating],
            r: float = wgs84.radiusearthkm) -> npt.NDArray[np.floating]:
        # Following https://en.wikipedia.org/wiki/Line–sphere_intersection
        # p0 is the satellite postion, with dimensions n_times, 3
        # p1 is the swath position, with dimensionsn_times, n_xtrack, 3

        # o is contains the same information as p_sat, buth with same
        # dimensions as p1
        o = np.ones(p.shape)
        o[:, :, 0] = (o[:, :, 0].T * p_sat[:, 0].T).T
        o[:, :, 1] = (o[:, :, 1].T * p_sat[:, 1].T).T
        o[:, :, 2] = (o[:, :, 2].T * p_sat[:, 2].T).T

        c = np.zeros(p.shape)

        # reshape o and p
        shp_ori = p.shape
        shp = (np.prod(shp_ori[:-1]), 3)
        o = (np.reshape(o, shp)).T
        p = (np.reshape(p, shp)).T
        c = (np.reshape(c, shp)).T

        # vector v from o to p and direction u
        v = p - o
        u = v / np.sqrt(np.sum(v**2, axis=0))

        u_dot_omc = np.sum(u * (o-c), axis=0)

        D = (u_dot_omc)**2 - np.sum((o-c) * (o-c), axis=0) + r**2
        D[D < 0.] = np.nan  # in these cases tthere are no intersects

        d = np.array([-u_dot_omc - np.sqrt(D), -u_dot_omc + np.sqrt(D)])

        p_int = o + np.min(d, axis=0) * u
        p_int = np.reshape(p_int.T, shp_ori)

        return p_int

    def compute_groundpoints(self,
                             satpos: dict,
                             yaw: float = 0.0,
                             pitch: float = 0.0,
                             roll: float = 0.0) -> dict:
        n_t = satpos['p'].shape[0]

        self.ground_points = {
            'seconds_from_start': satpos['seconds_from_start'],
            'thetas': self.thetas,
            'attitude': np.array([yaw, pitch, roll]),
        }

        # start with defining the unit vectors with satellite as
        # origin and pointing towards Earth center
        sat_x = (satpos['v'].T / np.sqrt(np.sum(satpos['v']**2, axis=-1))).T
        sat_z = (-satpos['p'].T / np.sqrt(np.sum(satpos['p']**2, axis=-1))).T
        sat_y = np.cross(sat_x, sat_z)
        sat_y = (sat_y.T / np.sqrt(np.sum(sat_y**2, axis=-1))).T

        # create look vectors
        # distance from the  Earth centre
        l_dist = np.sqrt(np.sum(satpos['p']**2, axis=-1))

        # first in the sat frame:
        p_s = np.zeros((n_t, self.thetas.shape[0], 3))
        for i, theta in enumerate(self.thetas):
            p_s[:, i, 1] = -l_dist * np.sin(np.deg2rad(theta))
            p_s[:, i, 2] = l_dist * np.cos(np.deg2rad(theta))

        # Apply in order roll-pitch-yaw. Change pitch sign to define
        # positive up.
        for i_t in range(n_t):
            A = Sensor.trans_yaw_pitch_roll(yaw, -pitch, roll, degrees=True)
            p_s[i_t, :, :] = np.matmul(p_s[i_t, :, :], A)

        # Translate to xyz axis
        p = np.zeros((n_t, self.thetas.shape[0], 3))
        for i in range(len(self.thetas)):
            p[:, i, :] = (p_s[:, i, 0] * sat_x.T + p_s[:, i, 1] * sat_y.T
                          + p_s[:, i, 2] * sat_z.T).T
            p[:, i, :] = satpos['p'][:, :] + p[:, i, :]

        # compute the intersection of the line between the satellite
        # and the points p and the spherical Earth
        self.ground_points['p'] = self.line_sphere_intersect(
            satpos['p'], p, r=wgs84.radiusearthkm)

        # compute geodetic info
        locs = EarthLocation.from_geocentric(
            self.ground_points['p'][:, :, 0],
            self.ground_points['p'][:, :, 1],
            self.ground_points['p'][:, :, 2],
            unit=u.km)
        self.ground_points['lat'] = locs.geodetic.lat.value
        self.ground_points['lon'] = locs.geodetic.lon.value
        self.ground_points['height'] = locs.geodetic.height.value

        # prepare spline interpolation of the points
        self.ground_points['time_scaled_axis'] = Sensor.scale_axis(
            self.ground_points['seconds_from_start'])
        self.ground_points['thetas_scaled_axis'] = Sensor.scale_axis(
            self.ground_points['thetas'])
        self.ground_points['p_scale_factor'] = 0.5 * (
            np.max(self.ground_points['p']) - np.min(self.ground_points['p']))

        self.ground_points['f_px'] = RectBivariateSpline(
            self.ground_points['time_scaled_axis'][0],
            self.ground_points['thetas_scaled_axis'][0],
            self.ground_points['p'][:, :, 0]
            / self.ground_points['p_scale_factor'])

        self.ground_points['f_py'] = RectBivariateSpline(
            self.ground_points['time_scaled_axis'][0],
            self.ground_points['thetas_scaled_axis'][0],
            self.ground_points['p'][:, :, 1]
            / self.ground_points['p_scale_factor'])

        self.ground_points['f_pz'] = RectBivariateSpline(
            self.ground_points['time_scaled_axis'][0],
            self.ground_points['thetas_scaled_axis'][0],
            self.ground_points['p'][:, :, 2]
            / self.ground_points['p_scale_factor'])

        self.ground_points['f_sat_px'] = interp1d(
            self.ground_points['time_scaled_axis'][0],
            satpos['p'][:, 0] / self.ground_points['p_scale_factor'],
            kind='cubic')

        self.ground_points['f_sat_py'] = interp1d(
            self.ground_points['time_scaled_axis'][0],
            satpos['p'][:, 1] / self.ground_points['p_scale_factor'],
            kind='cubic')

        self.ground_points['f_sat_pz'] = interp1d(
            self.ground_points['time_scaled_axis'][0],
            satpos['p'][:, 2] / self.ground_points['p_scale_factor'],
            kind='cubic')

        self.thetas_scaled = Sensor.scale_axis(self.thetas)

        return self.ground_points

    @staticmethod
    def get_theta_az(p0: npt.NDArray[np.floating],
                     p1: npt.NDArray[np.floating],
                     lat: npt.NDArray[np.floating],
                     lon: npt.NDArray[np.floating]) -> tuple[
                         npt.NDArray[np.floating], npt.NDArray[np.floating]]:
        # https://github.com/kg4sgp/look-angles/blob/master/sphere/\
        # lookangle_sphere.m
        if len(p0.shape) == 2:
            p0 = p0.reshape(p0.shape[0], 1, p0.shape[1])

        if len(lat.shape) == 1:
            lat = lat.reshape(lat.shape[0], 1)

        if len(lon.shape) == 1:
            lon = lon.reshape(lon.shape[0], 1)

        v = np.zeros(p0.shape)
        for i_x in range(v.shape[1]):
            v[:, i_x, :] = p1[:, :] - p0[:, i_x, :]

        # v = ( v.T /  np.sqrt(np.sum(v**2, axis=2)).T ).T # normalise to 1

        sin_lat = np.sin(np.deg2rad(lat))
        cos_lat = np.cos(np.deg2rad(lat))
        sin_lon = np.sin(np.deg2rad(lon))
        cos_lon = np.cos(np.deg2rad(lon))

        # Transform vector v to Topocentric Horizon
        rot_s = (sin_lat * cos_lon * v[:, :, 0]
                 + sin_lat * sin_lon * v[:, :, 1] - cos_lat * v[:, :, 2])
        rot_e = -1 * sin_lon * v[:, :, 0] + cos_lon * v[:, :, 1]
        rot_z = (cos_lat * cos_lon * v[:, :, 0]
                 + cos_lat * sin_lon * v[:, :, 1] + sin_lat * v[:, :, 2])

        r_range = np.sqrt(rot_s**2 + rot_e**2 + rot_z**2)

        # zenith angle
        theta = np.where(r_range < 1e-20, np.pi/2., np.arccos(rot_z/r_range))

        # azimuth
        # avoid divide by zero
        az = np.where(np.abs(rot_s) < 1e-20, np.pi/2.,
                      np.arctan(-1*(rot_e/rot_s)))
        az = np.where(rot_s > 0., az + np.pi, az)
        az = np.where(az < 0, az + 2*np.pi, az)

        return np.rad2deg(np.squeeze(theta)), np.rad2deg(np.squeeze(az))

    def get_groundpoints(self,
                         orbit_timestamps: npt.NDArray[np.datetime64],
                         dt_start: datetime.datetime,
                         dt_end: datetime.datetime,
                         dt_interval: float,
                         thetas: npt.NDArray[np.floating],
                         get_view_angles: bool = True,
                         get_solar_angles: bool = True) -> Geometry:
        # get ground points based on interpolation of computed points
        n_times = int((dt_end - dt_start).total_seconds() / dt_interval)

        # compute time axis as astropy object
        ap_ts = Time([dt_start + datetime.timedelta(seconds=float(sec))
                     for sec in np.arange(
                             0, n_times*dt_interval, dt_interval)])

        # determine the requested points in scaled axis units
        tsec = (
            (dt_start
             - orbit_timestamps[0].astype(datetime.datetime)).total_seconds()
            + np.arange(n_times) * dt_interval)

        tsec_scaled = Sensor.apply_scale_axis(
            tsec, self.ground_points['time_scaled_axis'])
        if np.max(np.abs(tsec_scaled)) > 1.:
            log.error('time out of range')
            exit(1)

        thetas_scaled = Sensor.apply_scale_axis(
            thetas, self.ground_points['thetas_scaled_axis'])

        # interpolate satellite position
        log.info('    satellite position')
        sat_p = np.zeros((tsec_scaled.shape[0], 3)) * np.nan
        sat_p[:, 0] = (self.ground_points['f_sat_px'](tsec_scaled)
                       * self.ground_points['p_scale_factor'])
        sat_p[:, 1] = (self.ground_points['f_sat_py'](tsec_scaled)
                       * self.ground_points['p_scale_factor'])
        sat_p[:, 2] = (self.ground_points['f_sat_pz'](tsec_scaled)
                       * self.ground_points['p_scale_factor'])

        # Interpolate points
        log.info('    satellite ground points')
        p = (np.zeros((tsec_scaled.shape[0], thetas_scaled.shape[0], 3))
             * np.nan)
        p[:, :, 0] = (self.ground_points['f_px'](tsec_scaled, thetas_scaled)
                      * self.ground_points['p_scale_factor'])
        p[:, :, 1] = (self.ground_points['f_py'](tsec_scaled, thetas_scaled)
                      * self.ground_points['p_scale_factor'])
        p[:, :, 2] = (self.ground_points['f_pz'](tsec_scaled, thetas_scaled)
                      * self.ground_points['p_scale_factor'])

        locs = EarthLocation.from_geocentric(
            p[:, :, 0], p[:, :, 1], p[:, :, 2], unit=u.km)

        lat = locs.geodetic.lat.value
        lon = locs.geodetic.lon.value

        # Get viewing angles
        log.info('    viewing angles')
        vza, vaa = Sensor.get_theta_az(p, sat_p, lat, lon)
        vaa = np.where(vaa > 180.0, vaa-360., vaa)

        # Get solar angles
        log.info('    solar angles')
        sun_p = get_sun(ap_ts).itrs.cartesian.xyz.value.T * u.AU.to(u.km)
        sza, saa = Sensor.get_theta_az(p, sun_p, lat, lon)
        saa = np.where(saa > 180.0, saa-360., saa)

        return Geometry(np.deg2rad(lat),
                        np.deg2rad(lon),
                        np.zeros(lat.shape),
                        np.deg2rad(sza),
                        np.deg2rad(saa),
                        np.deg2rad(vza),
                        np.deg2rad(vaa))
