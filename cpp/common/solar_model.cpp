// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "solar_model.h"

#include "constants.h"

namespace tango {

auto solarModel(const uint32_t tai_seconds,
                const double tai_second_fraction,
                Eigen::Vector3d& sun,
                Eigen::Quaterniond& q_j2000_ecef) -> void
{
    // UT1 time in modified Julian day 2000 (MJD2000) format
    double t_mjd {};
    // Fractional part of the number of UT1 days passed
    double day_fract_UT1 {};
    daysSinceJ2000(tai_seconds, tai_second_fraction, day_fract_UT1, t_mjd);

    // Compute solar ephemeris parameters. These are used to compute
    // the solar longitude and the nutation in longitude and
    // obliquity.
    double sun_mean_lon {};
    double sma {};
    solarEphemeris(t_mjd, sun_mean_lon, sma);
    const auto [moon_mean_lon, moon_asc_node] { lunarEphemeris(t_mjd) };

    // Compute nutation corrections for this day
    double nutation_lon {};
    double mean_obliquity {};
    double nutation_obl {};
    nutate(t_mjd,
           sun_mean_lon,
           sma,
           moon_mean_lon,
           moon_asc_node,
           nutation_lon,
           mean_obliquity,
           nutation_obl);

    // Greenwich hour angle in radians
    const double gha { getGHA2000(
      t_mjd, day_fract_UT1, nutation_lon, nutation_obl) };

    // Transform Greenwich hour angle to rotation Quaternion
    const Eigen::Quaterniond q_gha(
      std::cos(gha / 2.0), 0.0, 0.0, -std::sin(gha / 2.0));

    // Transform the sun and the satellite into ECR. The nutation and
    // such is already done. Only the Greenwich hour angle has to be
    // applied.
    sun = getSun2000(
      t_mjd, sun_mean_lon, sma, moon_mean_lon, nutation_lon, nutation_obl);
    sun = q_gha * sun;

    // J2000 to ECR rotation of the Sun
    Eigen::Quaterniond q_j2000_to_mod = j2000ToMOD(t_mjd);

    // Nutation quaternion
    Eigen::Quaterniond q_nutation =
      getNutationQuaternion(nutation_lon, nutation_obl, mean_obliquity);
    q_nutation = q_nutation.conjugate();

    // We can now compute the ECR transformation for the satellite
    // which is a combination of the Greenwich hour angle rotation,
    // nutation rotation, and the J2000 to MOD rotation.
    q_j2000_ecef = q_gha * q_nutation * q_j2000_to_mod;
}

[[nodiscard]] auto julianDayNumber(const int year,
                                   const int month,
                                   const int day) -> int
{
    const int a { (month - 14) / 12 };
    const int y { year + 4800 + a };
    return (1461 * y) / 4 + 367 * (month - 2 - 12 * a) / 12
           - (3 * ((y + 100) / 100)) / 4 + day - 32075;
}

[[nodiscard]] auto deltaTerrestrialUT1(const double y) -> double
{
    // https://eclipse.gsfc.nasa.gov/LEcat5/time.html
    if (y < 1961) {
        throw std::runtime_error { "years before 1961 not supported" };
    }
    if (y < 1986) {
        const double t { y - 1975.0 };
        const double t2 { t * t };
        const double t3 { t2 * t };
        return 45.45 + 1.067 * t - t2 / 260 - t3 / 718;
    }
    if (y < 2005) {
        const double t { y - 2000.0 };
        const double t2 { t * t };
        const double t3 { t2 * t };
        const double t4 { t3 * t };
        const double t5 { t4 * t };
        return 63.86 + 0.3345 * t - 0.060374 * t2 + 0.0017275 * t3
               + 0.000651814 * t4 + 0.00002373599 * t5;
    }
    if (y < 2050.0) {
        const double t { y - 2000.0 };
        const double t2 { t * t };
        return 62.92 + 0.32217 * t + 0.005589 * t2;
    }
    if (y < 2150) {
        const double d { ((y - 1820) / 100) };
        return -20 + 32 * d * d - 0.5628 * (2150 - y);
    }
    throw std::runtime_error { "years after 2150 not supported" };
}

auto daysSinceJ2000(const uint32_t tai_seconds,
                    const double tai_second_fraction,
                    double& day_fract_UT1,
                    double& t) -> void
{
    // Current Julian day number. The ICU time stamp is given relative
    // to 1 January 1958.
    const int days { static_cast<int>(tai_seconds / 86400) };
    const int julian_day_number { julianDayNumber(1958, 1, 1) + days };
    // Difference between terrestrial time (TT) and international
    // atomic time (TAI): TT = TAI + 32 s.
    constexpr double tai_tt_diff { 32.184 };
    // Number of years since 1958
    const double years { 1958 + tai_seconds / 86400 / 365.25 };
    // Difference between TAI and UT1
    const double delta_UT { deltaTerrestrialUT1(years) - tai_tt_diff };
    // Day fraction in UT1 time
    day_fract_UT1 =
      (tai_seconds % 86400 + tai_second_fraction - delta_UT) / 86400.0;
    const double julian_day { julian_day_number - 0.5 + day_fract_UT1 };
    t = julian_day - 2451545.0;
}

auto solarEphemeris(const double t, double& sun_mean_lon, double& sma) -> void
{
    // Sun mean longitude
    sun_mean_lon = fmod(280.46592 + 0.9856473516 * t, 360.0) * math::deg_to_rad;
    // Sun mean anomaly
    sma = fmod(357.52772 + 0.9856002831 * t, 360.0) * math::deg_to_rad;
}

[[nodiscard]] auto lunarEphemeris(const double t) -> std::array<double, 2>
{
    return { // Moon mean longitude
             fmod(218.31643 + 13.17639648 * t, 360.0) * math::deg_to_rad,
             // Ascending node of Moon's mean orbit
             fmod(125.04452 - 0.0529537648 * t, 360.0) * math::deg_to_rad
    };
}

auto nutate(const double t,
            const double sun_mean_lon,
            const double sma,
            const double moon_mean_lon,
            const double moon_asc_node,
            double& nutation_lon,
            double& mean_obliquity,
            double& nutation_obl) -> void
{
    // Nutation in Longitude (conversion from arcseconds to radians at the end)
    nutation_lon =
      (-17.1996 * std::sin(moon_asc_node)
       + 0.2062 * std::sin(2.0 * moon_asc_node)
       - 1.3187 * std::sin(2.0 * sun_mean_lon) + 0.1426 * std::sin(sma)
       - 0.2274 * std::sin(2.0 * moon_mean_lon))
      * math::deg_to_rad / 3600.0;
    //  Mean Obliquity of the Ecliptic
    mean_obliquity = (23.439291 - 3.560e-7 * t) * math::deg_to_rad;
    //  Nutation in Obliquity
    nutation_obl = mean_obliquity
                   + (9.2025 * std::cos(moon_asc_node)
                      + 0.5736 * std::cos(2.0 * sun_mean_lon))
                       * math::deg_to_rad / 3600.0;
}

auto getSun2000(const double t,
                const double sun_mean_lon,
                const double sma,
                const double moon_mean_lon,
                const double nutation_lon,
                const double nutation_obl) -> Eigen::Vector3d
{
    // Compute planet mean anomalies
    const double anomaly_venus { fmod(50.40828 + 1.60213022 * t, 360.0)
                                 * math::deg_to_rad };
    const double anomaly_mars { fmod(19.38816 + 0.52402078 * t, 360.0)
                                * math::deg_to_rad };
    const double anomaly_jupiter { fmod(20.35116 + 0.08309121 * t, 360.0)
                                   * math::deg_to_rad };
    // Compute geometric solar longitude
    const double geom_sun_lon {
        sun_mean_lon
        + ((6893.0 - 4.6543463e-4 * t) * std::sin(sma)
           + 72.0 * std::sin(2.0 * sma) - 7.0 * std::cos(sma - anomaly_jupiter)
           + 6.0 * std::sin(moon_mean_lon - sun_mean_lon)
           + 5.0
               * std::sin(4.0 * sma - 8.0 * anomaly_mars
                          + 3.0 * anomaly_jupiter)
           - 5.0 * std::cos(2.0 * sma - 2.0 * anomaly_venus)
           - 4.0 * std::sin(sma - anomaly_venus)
           + 4.0
               * std::cos(4.0 * sma - 8.0 * anomaly_mars
                          + 3.0 * anomaly_jupiter)
           + 3.0 * std::sin(2.0 * sma - 2.0 * anomaly_venus)
           - 3.0 * std::sin(anomaly_jupiter)
           - 3.0 * std::sin(2.0 * sma - 2.0 * anomaly_jupiter))
            * math::deg_to_rad / 3600.0
    };
    // Solar distance in astronomical units (AU)
    const double sun_dist { 1.00014 - 0.01671 * std::cos(sma)
                            - 0.00014 * std::cos(2.0 * sma) };
    // Constant of aberration
    const double abberation { 0.0056932 * math::deg_to_rad };
    // Apparent solar longitude including corrections for nutation in
    // longitude and velocity aberration
    const double app_sun_lon { geom_sun_lon + nutation_lon
                               - abberation / sun_dist };
    return { std::cos(app_sun_lon),
             std::sin(app_sun_lon) * std::cos(nutation_obl),
             std::sin(app_sun_lon) * std::sin(nutation_obl) };
}

[[nodiscard]] auto getGHA2000(const double t,
                              const double day_fract,
                              const double nutation_lon,
                              const double nutation_obl) -> double
{
    // Greenwich Mean Sidereal Time (radians)
    const double gmst { fmod(100.4606184 + 0.9856473663 * t + 2.908e-13 * t * t,
                             360.0)
                        * math::deg_to_rad };
    // Include apparent time correction and time-of-day
    return { gmst + nutation_lon * std::cos(nutation_obl)
             + day_fract * 360.0 * math::deg_to_rad };
}

auto j2000ToMOD(const double t) -> Eigen::Quaterniond
{
    // UT MJD2000 time in J2000 centuries
    const double tc { (t - 0.5) / 36525.0 };
    const double tc2 { tc * tc };
    const double tc3 { tc2 * tc };
    const double zeta0 { (2306.2181 * tc + 0.302 * tc2 + 0.018 * tc3)
                         * math::deg_to_rad / 3600.0 };
    const double thetap { (2004.3109 * tc - 0.4266 * tc2 - 0.04160 * tc3)
                          * math::deg_to_rad / 3600.0 };
    const double xip { (2306.2181 * tc + 1.095 * tc2 + 0.018 * tc3)
                       * math::deg_to_rad / 3600.0 };
    const Eigen::Matrix3d rot_mat {
        // 00
        { -std::sin(zeta0) * std::sin(xip)
            + std::cos(zeta0) * std::cos(xip) * std::cos(thetap),
          // 01
          -std::cos(zeta0) * std::sin(xip)
            - std::sin(zeta0) * std::cos(xip) * std::cos(thetap),
          // 02
          -std::cos(xip) * std::sin(thetap) },
        // 10
        { std::sin(zeta0) * std::cos(xip)
            + std::cos(zeta0) * std::sin(xip) * std::cos(thetap),
          // 11
          std::cos(zeta0) * std::cos(xip)
            - std::sin(zeta0) * std::sin(xip) * std::cos(thetap),
          // 12
          -std::sin(xip) * std::sin(thetap) },
        // 20
        { std::cos(zeta0) * std::sin(thetap),
          // 21
          -std::sin(zeta0) * std::sin(thetap),
          // 22
          std::cos(thetap) }
    };
    return Eigen::Quaterniond(rot_mat);
}

auto getNutationQuaternion(const double nutation_lon,
                           const double nutation_obl,
                           const double mean_obliquity) -> Eigen::Quaterniond
{
    // Nutation rotation matrix
    const Eigen::Matrix3d rot_mat {
        // 00
        { std::cos(nutation_lon),
          // 01
          std::sin(nutation_lon) * std::cos(nutation_obl),
          // 01
          std::sin(nutation_lon) * std::sin(nutation_obl) },
        // 10
        { -std::sin(nutation_lon) * std::cos(mean_obliquity),
          // 11
          std::cos(nutation_lon) * std::cos(nutation_obl)
              * std::cos(mean_obliquity)
            + std::sin(nutation_obl) * std::sin(mean_obliquity),
          // 12
          std::cos(nutation_lon) * std::sin(nutation_obl)
              * std::cos(mean_obliquity)
            - std::cos(nutation_obl) * std::sin(mean_obliquity) },
        // 20
        { -std::sin(nutation_lon) * std::sin(mean_obliquity),
          // 21
          std::cos(nutation_lon) * std::cos(nutation_obl)
              * std::sin(mean_obliquity)
            - std::sin(nutation_obl) * std::cos(mean_obliquity),
          // 22
          std::cos(nutation_lon) * std::sin(nutation_obl)
              * std::sin(mean_obliquity)
            + std::cos(nutation_obl) * std::cos(mean_obliquity) }
    };
    return Eigen::Quaterniond(rot_mat);
}

} // namespace tango
