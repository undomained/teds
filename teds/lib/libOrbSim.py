#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
   @author: Pepijn Veefkind   
            The Royal Netherlands Meteorological Institute (KNMI)
"""

import datetime
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import TEME, ITRS, CartesianDifferential, CartesianRepresentation, EarthLocation, SkyCoord, get_sun, AltAz
from sgp4.earth_gravity import wgs84
from sgp4.ext import jday
from sgp4.exporter import export_tle
from sgp4.model import Satrec, WGS84
from sgp4.api import SGP4_ERRORS


class Sensor:

    def trans_yaw_pitch_roll(self, y, p, r, degrees=True):
        # transformation matrix, see http://planning.cs.uiuc.edu/node102.html
        # ordering is: yaw, pitch, roll
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

        rot_mat = np.array([[cos_y*cos_p, cos_y*sin_p*sin_r - sin_y*cos_r, cos_y*sin_p*cos_r + sin_y*sin_r],
                            [sin_y*cos_p, sin_y*sin_p*sin_r + cos_y*cos_r, sin_y*sin_p*cos_r - cos_y*sin_r],
                            [-sin_p, cos_p*sin_r, cos_p*cos_r]])

        return rot_mat

    def scale_axis(self, x):
        # scale axis from -1 to 1
        offset = np.mean(x)
        scale_factor = 0.5 * (np.max(x) - np.min(x))
        return (x - offset) / scale_factor, offset, scale_factor

    def apply_scale_axis(self, x, scaled_axis):
        # apply scale factor and offset
        return (x - scaled_axis[1]) / scaled_axis[2]

    def apply_inv_scale_axis(self, x, scaled_axis):
        # apply scale factor and offset
        return x * scaled_axis[2] + scaled_axis[1]

    def dot(self, x, y):
        # function computes the dot product for the first dimesnion of a 2D matrix
        return(np.sum([x[i, :] * y[i, :] for i in range(x.shape[0])], axis=0))

    def R_wgs84_lat(self, lat, a=6378.1370, b=6356.7523142):
        # compute radius in km as a function of the geodetic latitude
        tan_beta = b / a * np.tan(np.deg2rad(lat))
        return np.sqrt((a**4 + b**4 * tan_beta**2) / (a**2 + b**2 * tan_beta**2))

    def line_sphere_intersect(self, p_sat, p, c=[0., 0., 0.], r=wgs84.radiusearthkm):
        # following https://en.wikipedia.org/wiki/Line–sphere_intersection
        # p0 is the satellite postion, with dimensions n_times, 3
        # p1 is the swath position, with dimensionsn_times, n_xtrack, 3

        # o is contains the same information as p_sat, buth with same dimensions as p1
        o = np.ones(p.shape)
        o[:, :, 0] = (o[:, :, 0].T * p_sat[:, 0].T).T
        o[:, :, 1] = (o[:, :, 1].T * p_sat[:, 1].T).T
        o[:, :, 2] = (o[:, :, 2].T * p_sat[:, 2].T).T

        c = np.zeros(p.shape) + c

        # reshape o and p
        shp_ori = p.shape
        shp = (np.prod(shp_ori[:-1]), 3)
        o = (np.reshape(o, shp)).T
        p = (np.reshape(p, shp)).T
        c = (np.reshape(c, shp)).T

        # vector v from o to p and direction u
        v = p - o
        u = v / np.sqrt(np.sum(v**2, axis=0))

        u_dot_omc = self.dot(u, (o-c))

        D = (u_dot_omc)**2 - self.dot((o-c), (o-c)) + r**2
        D[D < 0.] = np.nan  # in these cases tthere are no intersects

        d = np.array([-u_dot_omc - np.sqrt(D), -u_dot_omc + np.sqrt(D)])

        p_int = o + np.min(d, axis=0) * u
        p_int = np.reshape(p_int.T, shp_ori)

        return(p_int)

    def compute_groundpoints(self, satpos, yaw=0., pitch=0., roll=0.):

        n_t = satpos['p'].shape[0]

        self.ground_points = {}
        self.ground_points['datetime'] = satpos['datetime']
        self.ground_points['seconds_from_start'] = satpos['seconds_from_start']
        self.ground_points['thetas'] = self.thetas
        self.ground_points['attitude'] = np.array([yaw, pitch, roll])

        # start with defining the unit vectors with satellite as origin and pointing towards Earth center
        sat_x = (satpos['v'].T / np.sqrt(np.sum(satpos['v']**2, axis=-1))).T
        sat_z = (-satpos['p'].T / np.sqrt(np.sum(satpos['p']**2, axis=-1))).T
        sat_y = np.cross(sat_x, sat_z)
        sat_y = (sat_y.T / np.sqrt(np.sum(sat_y**2, axis=-1))).T

        # create look vectors
        l = np.sqrt(np.sum(satpos['p']**2, axis=-1))  # distance from the  Earth centre

        # first in the sat frame:
        p_s = np.zeros((n_t, self.thetas.shape[0], 3))
        for i, theta in enumerate(self.thetas):
            p_s[:, i, 1] = -l * np.sin(np.deg2rad(theta))
            p_s[:, i, 2] = l * np.cos(np.deg2rad(theta))

        # apply in order roll-pitch-yaw
        # change pitch sign to define positive up
        for i_t in range(n_t):
            A = self.trans_yaw_pitch_roll(yaw, -pitch, roll, degrees=True)
            p_s[i_t, :, :] = np.matmul(p_s[i_t, :, :], A)

        # translate to xyz axis
        p = np.zeros((n_t, self.thetas.shape[0], 3))
        for i, theta in enumerate(self.thetas):
            p[:, i, :] = (p_s[:, i, 0] * sat_x.T + p_s[:, i, 1] * sat_y.T + p_s[:, i, 2] * sat_z.T).T
            p[:, i, :] = satpos['p'][:, :] + p[:, i, :]

        # compute the intersection of the line between the satellite and the points p and the spherical Earth
        self.ground_points['p'] = self.line_sphere_intersect(satpos['p'], p, c=[0., 0., 0.], r=wgs84.radiusearthkm)

        # compute geodetic info
        locs = EarthLocation.from_geocentric(
            self.ground_points['p'][:, :, 0], self.ground_points['p'][:, :, 1], self.ground_points['p'][:, :, 2], unit=u.km)
        self.ground_points['lat'] = locs.geodetic.lat.value
        self.ground_points['lon'] = locs.geodetic.lon.value
        self.ground_points['height'] = locs.geodetic.height.value

        # prepare spline interpolation of the points
        self.ground_points['time_scaled_axis'] = self.scale_axis(self.ground_points['seconds_from_start'])
        self.ground_points['thetas_scaled_axis'] = self.scale_axis(self.ground_points['thetas'])
        self.ground_points['p_scale_factor'] = 0.5 * (np.max(self.ground_points['p']) - np.min(self.ground_points['p']))

        self.ground_points['f_px'] = RectBivariateSpline(self.ground_points['time_scaled_axis'][0],
                                                         self.ground_points['thetas_scaled_axis'][0],
                                                         self.ground_points['p'][:, :, 0] / self.ground_points['p_scale_factor'])

        self.ground_points['f_py'] = RectBivariateSpline(self.ground_points['time_scaled_axis'][0],
                                                         self.ground_points['thetas_scaled_axis'][0],
                                                         self.ground_points['p'][:, :, 1] / self.ground_points['p_scale_factor'])

        self.ground_points['f_pz'] = RectBivariateSpline(self.ground_points['time_scaled_axis'][0],
                                                         self.ground_points['thetas_scaled_axis'][0],
                                                         self.ground_points['p'][:, :, 2] / self.ground_points['p_scale_factor'])

        self.ground_points['f_sat_px'] = interp1d(self.ground_points['time_scaled_axis'][0],
                                                  satpos['p'][:, 0] / self.ground_points['p_scale_factor'], kind='cubic')

        self.ground_points['f_sat_py'] = interp1d(self.ground_points['time_scaled_axis'][0],
                                                  satpos['p'][:, 1] / self.ground_points['p_scale_factor'], kind='cubic')

        self.ground_points['f_sat_pz'] = interp1d(self.ground_points['time_scaled_axis'][0],
                                                  satpos['p'][:, 2] / self.ground_points['p_scale_factor'], kind='cubic')

        self.thetas_scaled = self.scale_axis(self.thetas)

        return self.ground_points

    def get_solar_angles(self, ap_ts, locs):
        # This function uses astropy and includes refraction in the atmosphere.
        # However, it is much slower compared to get_theta_az.

        suns = get_sun(ap_ts)

        sza = np.zeros(locs.shape) * np.nan
        saa = np.zeros(locs.shape) * np.nan
        for i, sun in enumerate(suns):
            altaz = sun.transform_to(AltAz(location=locs[i, :]))

            saa[i, :] = altaz.az.value
            sza[i, :] = 90. - altaz.alt.value

        return sza, saa

    def get_viewing_angles(self, ap_ts, locs, sat_locs):
        # This function uses astropy and includes refraction in the atmosphere.
        # However, it is much slower compared to get_theta_az.

        vza = np.zeros(locs.shape) * np.nan
        vaa = np.zeros(locs.shape) * np.nan

        for i, sat_loc in enumerate(sat_locs):

            sc_sat_loc = SkyCoord(sat_loc.get_gcrs(ap_ts[i]))

            altaz = sc_sat_loc.transform_to(AltAz(location=locs[i, :]))

            vaa[i, :] = altaz.az.value
            vza[i, :] = 90. - altaz.alt.value

        return vza, vaa

    def get_theta_az(self, p0, p1, lat, lon):
        # using method from https://github.com/kg4sgp/look-angles/blob/master/sphere/lookangle_sphere.m
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
        rot_s = sin_lat * cos_lon * v[:, :, 0] + sin_lat * sin_lon * v[:, :, 1] - cos_lat * v[:, :, 2]
        rot_e = -1 * sin_lon * v[:, :, 0] + cos_lon * v[:, :, 1]
        rot_z = cos_lat * cos_lon * v[:, :, 0] + cos_lat * sin_lon * v[:, :, 1] + sin_lat * v[:, :, 2]

        r_range = np.sqrt(rot_s**2 + rot_e**2 + rot_z**2)

        # zenith angle
        theta = np.where(r_range < 1e-20, np.pi/2., np.arccos(rot_z/r_range))

        # azimuth
        az = np.where(np.abs(rot_s) < 1e-20, np.pi/2., np.arctan(-1*(rot_e/rot_s)))  # avoid divide by zero
        az = np.where(rot_s > 0., az + np.pi, az)
        az = np.where(az < 0, az + 2*np.pi, az)

        return np.rad2deg(np.squeeze(theta)), np.rad2deg(np.squeeze(az))

    def get_groundpoints(self, dt_start, dt_end, dt_interval, thetas, get_view_angles=True, get_solar_angles=True):
        # get ground points based on interpolation of computed points
        n_times = int((dt_end - dt_start).total_seconds() / dt_interval) + 1

        # compute time axis as astropy object
        ap_ts = Time([dt_start + datetime.timedelta(seconds=sec)
                     for sec in np.arange(0, n_times*dt_interval, dt_interval)])

        # determine the requested points in scaled axis units
        tsec = (dt_start - self.ground_points['datetime'][0]).total_seconds() + np.arange(n_times) * dt_interval

        tsec_scaled = self.apply_scale_axis(tsec, self.ground_points['time_scaled_axis'])
        if np.max(np.abs(tsec_scaled)) > 1.:
            print('ERROR: time out of range')
            return None

        thetas_scaled = self.apply_scale_axis(thetas, self.ground_points['thetas_scaled_axis'])
        if np.max(np.abs(thetas_scaled)) > 1.:
            print('ERROR: theta  out of range')
            return None

        # interpolate satellite position
        print('    satellite position')
        sat_p = np.zeros((tsec_scaled.shape[0], 3)) * np.nan
        sat_p[:, 0] = self.ground_points['f_sat_px'](tsec_scaled) * self.ground_points['p_scale_factor']
        sat_p[:, 1] = self.ground_points['f_sat_py'](tsec_scaled) * self.ground_points['p_scale_factor']
        sat_p[:, 2] = self.ground_points['f_sat_pz'](tsec_scaled) * self.ground_points['p_scale_factor']

        sat_locs = EarthLocation.from_geocentric(sat_p[:, 0], sat_p[:, 1], sat_p[:, 2], unit=u.km)
        sat_lat = sat_locs.geodetic.lat.value
        sat_lon = sat_locs.geodetic.lon.value
        sat_height = sat_locs.geodetic.height.value

        # interpolate points
        print('    satellite ground points')
        p = np.zeros((tsec_scaled.shape[0], thetas_scaled.shape[0], 3)) * np.nan
        p[:, :, 0] = self.ground_points['f_px'](tsec_scaled, thetas_scaled) * self.ground_points['p_scale_factor']
        p[:, :, 1] = self.ground_points['f_py'](tsec_scaled, thetas_scaled) * self.ground_points['p_scale_factor']
        p[:, :, 2] = self.ground_points['f_pz'](tsec_scaled, thetas_scaled) * self.ground_points['p_scale_factor']

        locs = EarthLocation.from_geocentric(p[:, :, 0], p[:, :, 1], p[:, :, 2], unit=u.km)

        lat = locs.geodetic.lat.value
        lon = locs.geodetic.lon.value
        height = locs.geodetic.height.value

        # get viewing angles
        if get_view_angles:
            print('    viewing angles')
            vza, vaa = self.get_theta_az(p, sat_p, lat, lon)
            vaa = np.where(vaa > 180.0, vaa-360., vaa)  # scale between -180 and 180
        else:
            vza = None
            vaa = None

        # get solar angles
        if get_solar_angles:
            print('    solar angles')
            sun_p = get_sun(ap_ts).itrs.cartesian.xyz.value.T * u.AU.in_units('km')
            sza, saa = self.get_theta_az(p, sun_p, lat, lon)
            saa = np.where(saa > 180.0, saa-360., saa)  # scale between -180 and 180
        else:
            sza = None
            saa = None

        return {'p': p, 'lat': lat, 'lon': lon, 'height': height,
                'sat_p': sat_p, 'sat_lat': sat_lat, 'sat_lon': sat_lon, 'sat_height': sat_height,
                'start_time': dt_start, 'end_time': dt_end, 'seconds_from_epoch': tsec,
                'epoch': self.ground_points['datetime'][0],
                'vza': vza, 'vaa': vaa, 'sza': sza, 'saa': saa,
                'thetas': thetas}

    def __init__(self, theta_start, theta_end, d_theta):

        self.thetas = np.arange(theta_start, theta_end + d_theta, d_theta)
        self.n_thetas = self.thetas.shape[0]

        return


class Satellite:

    # Constants
    grav_model = WGS84
    mu = wgs84.mu                        # km^3  s^-2
    earth_radius = wgs84.radiusearthkm   # km
    j2 = wgs84.j2

    opsmode = 'i'

    def doy(self, dt):
        # return day of year as a float
        dt_doy = dt-datetime.datetime(dt.year, 1, 1, 0, 0, 0)
        return dt_doy.days + (dt_doy.seconds + 1e-6*dt_doy.microseconds) / (60.*60.*24)

    def sun_synchronous(self, h, ltan, eccentricity=1e-9):
        # height in km
        # ltan in hours

        # semi-major axis
        sma = h + self.earth_radius

        # mean motion
        mean_motion = np.sqrt(self.mu / sma**3)               # radians / second

        inclination = np.arccos(-1 * (sma / 12352.0)**(7./2.))  # radians

        # RAAN
        epoch = Time(self.epoch)  # copy to astropy object
        raan = 2*np.pi * (epoch.sidereal_time("mean", 0).value / 24. +
                          ((0.5 - epoch.jd2) + ltan / 24.))  # radians
        raan = raan % (2*np.pi)

        return sma, inclination, mean_motion, raan

    def compute_position(self, dt_start, dt_end, dt_interval=10.0):
        # compute the number of times
        n = int((dt_end - dt_start).total_seconds() / dt_interval) + 1

        dts = [dt_start + datetime.timedelta(seconds=sec) for sec in np.arange(0, n*dt_interval, dt_interval)]

        ap_ts = Time(dts)

        pos = {}
        pos['datetime'] = np.array(dts)
        pos['seconds_from_start'] = np.arange(0, n*dt_interval, dt_interval)
        pos['lat'] = np.zeros(n) * np.nan
        pos['lon'] = np.zeros(n) * np.nan
        pos['height'] = np.zeros(n) * np.nan
        pos['p'] = np.zeros((n, 3)) * np.nan
        pos['v'] = np.zeros((n, 3)) * np.nan

        for i, t in enumerate(ap_ts):
            error_code, teme_p, teme_v = self.satrec.sgp4(t.jd1, t.jd2)  # in km and km/s
            if error_code != 0:
                raise RuntimeError(SGP4_ERRORS[error_code])

            teme_p = CartesianRepresentation(teme_p*u.km)
            teme_v = CartesianDifferential(teme_v*u.km/u.s)
            teme = TEME(teme_p.with_differentials(teme_v), obstime=t)

            itrs = teme.transform_to(ITRS(obstime=t))
            loc = itrs.earth_location

            pos['lat'][i] = loc.geodetic.lat.value
            pos['lon'][i] = loc.geodetic.lon.value
            pos['height'][i] = loc.geodetic.height.value

            pos['p'][i, :] = [itrs.x.value, itrs.y.value, itrs.z.value]
            pos['v'][i, :] = [itrs.v_x.value, itrs.v_y.value, itrs.v_z.value]

        return pos

    def __init__(self, orbdef, satnum=99999, nddot=0.0, intldesg='12345A', classification='U', revnum=0, elnum=1):

        # init internal satrec object
        self.satrec = Satrec()

        # define the epoch year and doy
        dt = orbdef['epoch']
        self.satrec.jdsatepoch = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

        # for convenience, also store the epoch as datetime object
        self.epoch = orbdef['epoch']

        # define name and number
        if 'satnum' in orbdef.keys():
            self.satrec.satnum_str = orbdef['satnum']
        else:
            self.satrec.satnum_str = str(satnum)

        # define nddot
        if 'nddot' in orbdef.keys():
            self.satrec.nddot = orbdef['nddot']
        else:
            self.satrec.nddot = nddot

        # compute the inclination from the height
        #a, inclination, mean_motion = self.sun_synchronous(orbdef['sat_height'])
        sma, inclination, mean_motion, raan = self.sun_synchronous(orbdef['sat_height'], orbdef['ltan'])

        self.satrec.inclo = inclination
        self.satrec.no_kozai = mean_motion * 60.  # radians / minute
        self.satrec.nodeo = raan

        # copy orbital parameters
        self.satrec.ecco = orbdef['eccentricity']
        self.satrec.argpo = orbdef['arg_perigee']
        self.satrec.mo = orbdef['mean_anomaly']
        self.satrec.ndot = orbdef['mean_motion_dot']
        self.satrec.bstar = orbdef['drag_coeff']

        self.satrec.sgp4init(self.grav_model, self.opsmode, self.satrec.satnum,
                             self.satrec.jdsatepoch-2433281.5, self.satrec.bstar,
                             self.satrec.ndot, self.satrec.nddot, self.satrec.ecco,
                             self.satrec.argpo, self.satrec.inclo, self.satrec.mo,
                             self.satrec.no_kozai, self.satrec.nodeo)

        # load other TLE info
        if 'classification' in orbdef.keys():
            self.satrec.classification = orbdef['classification']
        else:
            self.satrec.classification = classification

        if 'intldesg' in orbdef.keys():
            self.satrec.intldesg = orbdef['intldesg']
        else:
            self.satrec.intldesg = intldesg

        if 'revnum' in orbdef.keys():
            self.satrec.revnum = orbdef['revnum']
        else:
            self.satrec.revnum = revnum

        self.satrec.ephtype = "0"  # always 0 in distributed ephemeris

        if 'elnum' in orbdef.keys():
            self.satrec.elnum = orbdef['elnum']
        else:
            self.satrec.elnum = elnum

        # compute the tle
        self.tle = export_tle(self.satrec)

        return


