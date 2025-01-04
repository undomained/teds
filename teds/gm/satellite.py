# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""SGP4 wrapper to generate orbit positions."""
from astropy import units as u
from astropy.coordinates import CartesianDifferential
from astropy.coordinates import CartesianRepresentation
from astropy.coordinates import ITRS
from astropy.coordinates import TEME
from astropy.time import Time
from sgp4.api import SGP4_ERRORS
from sgp4.earth_gravity import wgs84
from sgp4.exporter import export_tle
from sgp4.ext import jday
from sgp4.model import Satrec
from sgp4.model import WGS84
import numpy as np
import numpy.typing as npt


class Satellite:
    grav_model = WGS84
    mu = wgs84.mu                        # km^3  s^-2
    earth_radius = wgs84.radiusearthkm   # km
    j2 = wgs84.j2
    opsmode = 'i'

    def __init__(self, orbdef: dict) -> None:
        # Init internal satrec object
        self.satrec = Satrec()

        # Define the epoch year and doy
        dt = orbdef['epoch']
        self.satrec.jdsatepoch = jday(dt.year,
                                      dt.month,
                                      dt.day,
                                      dt.hour,
                                      dt.minute,
                                      dt.second)

        # For convenience, also store the epoch as datetime object
        self.epoch = orbdef['epoch']

        # Define name and number
        if 'satnum' in orbdef:
            self.satrec.satnum_str = str(orbdef['satnum'])
        else:
            self.satrec.satnum_str = '99999'

        # Define nddot
        if 'nddot' in orbdef:
            self.satrec.nddot = orbdef['nddot']
        else:
            self.satrec.nddot = 0.0

        # Compute the inclination from the height
        sma, inclination, mean_motion, raan = self.sun_synchronous(
            orbdef['sat_height'], orbdef['ltan'])

        self.satrec.inclo = inclination
        self.satrec.no_kozai = mean_motion * 60.  # radians / minute
        self.satrec.nodeo = raan

        # Copy orbital parameters
        self.satrec.ecco = orbdef['eccentricity']
        self.satrec.argpo = orbdef['arg_perigee']
        self.satrec.mo = orbdef['mean_anomaly']
        self.satrec.ndot = orbdef['mean_motion_dot']
        self.satrec.bstar = orbdef['drag_coeff']

        self.satrec.sgp4init(self.grav_model,
                             self.opsmode,
                             self.satrec.satnum,
                             self.satrec.jdsatepoch-2433281.5,
                             self.satrec.bstar,
                             self.satrec.ndot,
                             self.satrec.nddot,
                             self.satrec.ecco,
                             self.satrec.argpo,
                             self.satrec.inclo,
                             self.satrec.mo,
                             self.satrec.no_kozai,
                             self.satrec.nodeo)

        # Load other TLE info
        if 'classification' in orbdef:
            self.satrec.classification = orbdef['classification']
        else:
            self.satrec.classification = 'U'

        if 'intldesg' in orbdef:
            self.satrec.intldesg = orbdef['intldesg']
        else:
            self.satrec.intldesg = '12345A'

        if 'revnum' in orbdef:
            self.satrec.revnum = orbdef['revnum']
        else:
            self.satrec.revnum = 0

        self.satrec.ephtype = "0"  # always in distributed ephemeris

        if 'elnum' in orbdef:
            self.satrec.elnum = orbdef['elnum']
        else:
            self.satrec.elnum = 1

        # Compute the TLE
        self.tle = export_tle(self.satrec)

    def sun_synchronous(
            self,
            h: float,
            ltan: float,
            eccentricity: float = 1e-9) -> tuple[float, float, float, float]:
        # height in km
        # ltan in hours

        # semi-major axis
        sma = h + self.earth_radius

        # mean motion
        mean_motion = np.sqrt(self.mu / sma**3)  # radians / second

        inclination = np.arccos(-1 * (sma / 12352.0)**(7./2.))  # radians

        # RAAN
        epoch = Time(self.epoch)  # copy to astropy object
        raan = 2*np.pi * (epoch.sidereal_time("mean", 0).value / 24. +
                          ((0.5 - epoch.jd2) + ltan / 24.))  # radians
        raan = raan % (2*np.pi)

        return sma, inclination, mean_motion, raan

    def compute_positions(self,
                          timestamps: npt.NDArray[np.datetime64]) -> dict:
        pos = {
            'seconds_from_start': np.array(
                1e-6 * (timestamps - timestamps[0]), dtype=np.float64),
            'lat': np.zeros(len(timestamps)),
            'lon': np.zeros(len(timestamps)),
            'height': np.zeros(len(timestamps)),
            'p': np.zeros((len(timestamps), 3)),
            'p_j2000': np.zeros((len(timestamps), 3)),
            'v': np.zeros((len(timestamps), 3)),
        }
        for i, t in enumerate(Time(timestamps)):
            # in km and km/s
            error_code, teme_p, teme_v = self.satrec.sgp4(t.jd1, t.jd2)
            if error_code != 0:
                raise RuntimeError(SGP4_ERRORS[error_code])
            teme_p = CartesianRepresentation(teme_p * u.km)
            teme_v = CartesianDifferential(teme_v * u.km / u.s)
            teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
            itrs = teme.transform_to(ITRS(obstime=t))
            loc = itrs.earth_location
            pos['lat'][i] = loc.geodetic.lat.value
            pos['lon'][i] = loc.geodetic.lon.value
            pos['height'][i] = loc.geodetic.height.value
            pos['p'][i, :] = [itrs.x.value, itrs.y.value, itrs.z.value]
            pos['v'][i, :] = [itrs.v_x.value, itrs.v_y.value, itrs.v_z.value]
        return pos
