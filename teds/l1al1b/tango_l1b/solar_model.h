// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// The solar model. Much of it comes from the Astronomical Almanac
// (1983 or 1984).
//
// Calculate the sun position and the Greenwich hour angle. Both need
// a complicated astronomical model from the Astronomical Almanac. To
// keep the perspective as consistent as possible, we will calculate
// the sun position in J2000 inside the routine because the satellite
// position is also in J2000. The solar model works for Earth only.
//
// The Sun vector is computed in geocentric inertial (equatorial)
// coodinates. The accuracy of the Sun vector is approximately 0.1
// arcminute.

#pragma once

#include "constants.h"
#include "quaternion.h"

#include <array>
#include <cstdint>

namespace tango {

// The solar takes only one type of argument: the number of seconds
// since 01/01/1958 in the internation atomic time (TAI) format and
// the TAI second fractional part. The output is the sun vector in the
// ECR coordinates and the J200-to-ECR quaternion.
auto solarModel(const uint32_t tai_seconds,
                const double tai_second_fraction,
                std::array<double, dims::vec>& sun,
                Quaternion& q_j2000_ecef) -> void;

extern "C" auto c_solarModel(const uint32_t tai_time,
                             const double tai_time_subsec,
                             double* p_q_j2000_ecef) -> void;

// Number of days since noon January 1, 4713 BC of the Julian
// calendar. Arguments form a Gregorian calendar date.
[[nodiscard]] auto julianDayNumber(const int year,
                                   const int month,
                                   const int day) -> int;

// Return the time difference delta = TT - UT1, where TT is
// terrestrial time (y) and UT1 is universal time.
[[nodiscard]] auto deltaTerrestrialUT1(const double y) -> double;

// Return the day fraction in the UT1 time scale and the UT1 timestamp
// in the modified Julian day 2000 format (MJD2000).
auto daysSinceJ2000(const uint32_t tai_seconds,
                    const double tai_second_fraction,
                    double& day_fract_UT1,
                    double& t) -> void;

// Return the Sun mean longitude and anomaly. The formulas are in
// degrees and the modulo over 360 is calculated in advance to stay
// consistent with common implementation. Changing the order of modulo
// and unit conversion, or just execute sine or cosine without modulo
// over full circles leads to different numeric errors.
auto solarEphemeris(const double t, double& sun_mean_lon, double& sma) -> void;

// The Moon mean longitude and the ascending node of Moon's mean orbit
[[nodiscard]] auto lunarEphemeris(const double t) -> std::array<double, 2>;

// Compute the nutation in longitude and the obliquity of the ecliptic
// corrected for nutation. These parameters are used to compute the
// apparent time correction to the Greenwich Hour Angle and for the
// calculation of the geocentric Sun vector. Terms are included to 0.1
// arcsecond. Arguments are UT1 time in MJD2000, sun mean longitude,
// sun mean anomaly, moon mean longitude, and the ascending node of
// Moon's mean orbit.
auto nutate(const double t,
            const double sun_mean_lon,
            const double sma,
            const double moon_mean_lon,
            const double moon_asc_node,
            double& nutation_lon,
            double& mean_obliquity,
            double& nutation_obl) -> void;

// Compute the sun vector in geocentric inertial
// coordinates. Arguments are UT1 time in MJD2000, sun mean longitude,
// sun mean anomaly, moon mean longitude, nutation in Longitude, and
// Nutation in obliquity.
auto getSun2000(const double t,
                const double sun_mean_lon,
                const double sma,
                const double moon_mean_lon,
                const double nutation_lon,
                const double nutation_obl,
                std::array<double, dims::vec>& sun) -> void;

// Greenwich hour angle in radians
[[nodiscard]] auto getGHA2000(const double t,
                              const double day_fract,
                              const double nutation_lon,
                              const double nutation_obl) -> double;

// Get J2000 to mean-of-day (precession) transformation
auto j2000ToMOD(const double t, Quaternion& q_j2000_to_mod) -> void;

auto getNutationQuaternion(const double nutation_lon,
                           const double nutation_obl,
                           const double mean_obliquity,
                           Quaternion& q_nutation) -> void;

} // namespace tango
