// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "geometry.h"

#include "algorithm.h"
#include "solar_model.h"

#include <cmath>
#include <stdexcept>

namespace tango {

auto moduloAndConvertAngles(std::vector<double>& angles) -> void
{
    for (auto& angle : angles) {
        // This does nothing for latitudes
        if (angle > std::numbers::pi) {
            angle -= 2 * std::numbers::pi;
        }
        if (angle < -std::numbers::pi) {
            angle += 2 * std::numbers::pi;
        }
        angle /= static_cast<double>(math::deg_to_rad);
    }
}

auto interpolateQuaternion(const std::vector<double>& att_times,
                           const std::vector<Quaternion>& att_quaternions,
                           const double image_time,
                           Quaternion& quat) -> void
{
    // Find the first attitude time higher than the current time stamp
    const auto first_t_higher { std::ranges::find_if(
      att_times, [image_time](const double t) { return t > image_time; }) };
    // Index of the first attitude time lower than current time stamp
    const auto idx { std::distance(att_times.begin(), first_t_higher) - 1 };
    const double idx_delta { (image_time - att_times[idx])
                             / (att_times[idx + 1] - att_times[idx]) };
    quat = att_quaternions[idx];
    quat.slerp(att_quaternions[idx + 1], idx_delta);
}

auto intersectEllipsoid(std::array<double, dims::vec> los,
                        std::array<double, dims::vec> sat,
                        std::array<double, dims::vec>& pos) -> void
{
    los[0] /= earth::a;
    los[1] /= earth::a;
    los[2] /= earth::b;
    sat[0] /= earth::a;
    sat[1] /= earth::a;
    sat[2] /= earth::b;
    const double ps { dotProduct(los, sat) };
    const double pp { dotProduct(los, los) };
    const double ss { dotProduct(sat, sat) };
    const double discriminant { ps * ps - pp * (ss - 1.0) };
    if (discriminant < 0.0) {
        throw std::runtime_error { "no intersection with the ellipsoid, the "
                                   "satellite is looking past Earth" };
    }
    const double d { (-ps - std::sqrt(discriminant)) / pp };
    for (int i {}; i < dims::vec; ++i) {
        pos[i] = sat[i] + d * los[i];
    }
    pos[0] *= earth::a;
    pos[1] *= earth::a;
    pos[2] *= earth::b;
}

auto intersectTerrain(DEM& dem,
                      const std::array<double, dims::vec>& los,
                      std::array<double, dims::vec>& pos,
                      double& lat,
                      double& lon,
                      double& height) -> void
{
    cart2geo(pos, lat, lon, height);
    // Normal vector n at the target point
    const std::array<double, dims::vec> normal { std::cos(lat) * std::cos(lon),
                                                 std::cos(lat) * std::sin(lon),
                                                 std::sin(lat) };
    const double cosn { -dotProduct(los, normal) };
    double sign { -1.0 };
    // Highest point on Earth
    constexpr double max_jump { 10000.0 };
    double jump_size { max_jump / cosn };
    // Stop the algorithm if the altitude changes by less than this
    // distance in m.
    constexpr double target_accuracy { 1.0 };
    // Ideally the algorithm stops when the target accuracy is
    // reached. The maximum number of iterations is here to prevent an
    // infinite loop if something goes wrong. At this iteration level
    // the error is on the order of 1 cm.
    constexpr int max_iter { 20 };
    for (int iter {}; iter < max_iter; ++iter) {
        for (int i {}; i < dims::vec; ++i) {
            pos[i] += sign * jump_size * los[i];
        }
        const double height_prev { height };
        cart2geo(pos, lat, lon, height);
        if (std::abs(height_prev - height) < target_accuracy) {
            break;
        }
        // The next jump direction is determined by whether the
        // current point is below or above the height given by the DEM
        // file.
        sign = (dem.getHeight(lat, lon) < height ? 1.0 : -1.0);
        jump_size /= 2;
    }
}

auto cart2geo(const std::array<double, dims::vec>& xyz,
              double& lat,
              double& lon,
              double& height) -> void
{
    // Longitude is easy
    lon = std::atan2(xyz[1], xyz[0]);
    // Latitude and height are computed using the algorithm of
    // Fukushima "Transformation from Cartesian to Geodetic
    // Coordinates Accelerated by Halley's Method", J. Geodesy 79, 689
    // (2006).
    static constexpr double ec { std::sqrt(1 - earth::e2) };
    static constexpr double c { earth::a * earth::e2 };
    const double p { std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]) };
    const double zc { ec * std::abs(xyz[2]) };
    double Sn { std::abs(xyz[2]) };
    double Cn { ec * p };
    const double Sn_prev { Sn };
    // This algorithm is only run for one iteration
    const double SnCn { (Sn * Sn + Cn * Cn) };
    const double An { std::sqrt(SnCn) };
    const double cSnCn { c * Sn * Cn };
    const double Bn { 1.5 * cSnCn * ((p * Sn - zc * Cn) * An - cSnCn) };
    const double An3 { SnCn * An };
    const double Sn3 { Sn * Sn * Sn };
    const double Cn3 { Cn * Cn * Cn };
    Sn = (zc * An3 + c * Sn3) * An3 - Bn * Sn;
    Cn = (p * An3 - c * Cn3) * An3 - Bn * Cn;
    // End of iteration one
    const double Cc { ec * Cn };
    lat = std::atan(Sn / Cc);
    height = (p * Cc + Sn_prev * Sn
              - earth::a * std::sqrt(ec * ec * Sn * Sn + Cc * Cc))
             / (std::sqrt(Sn * Sn + Cc * Cc));
    if (xyz[2] < 0.0) {
        lat = -lat;
    }
}

auto solarAndViewingGeometry(const std::array<double, dims::vec>& sun,
                             const std::array<double, dims::vec>& los,
                             const double lat,
                             const double lon,
                             double& sza,
                             double& saa,
                             double& vza,
                             double& vaa) -> void
{
    // Normal vector of the target point
    const std::array<double, dims::vec> zenith { std::cos(lat) * std::cos(lon),
                                                 std::cos(lat) * std::sin(lon),
                                                 std::sin(lat) };
    // Normal vector in the east direction
    const std::array<double, dims::vec> east { -std::sin(lon),
                                               std::cos(lon),
                                               0 };
    // Normal vector in the north direction
    std::array<double, dims::vec> north {};
    crossProduct(zenith, east, north);
    // Solar geometry
    sza = std::acos(dotProduct(zenith, sun));
    saa = std::atan2(dotProduct(east, sun), dotProduct(north, sun));
    // Viewing geometry
    vza = std::acos(-dotProduct(zenith, los));
    vaa = std::atan2(-dotProduct(east, los), -dotProduct(north, los));
}

auto geolocate(const std::string& dem_filename,
               const std::vector<double>& los,
               const std::vector<uint32_t>& tai_seconds,
               const std::vector<double>& tai_subsec,
               const std::vector<double>& orb_pos_j2000,
               const std::vector<Quaternion>& att_quat_sc_j2000,
               Geometry& geo) -> void
{
    DEM dem { dem_filename };
    const int n_alt { static_cast<int>(tai_seconds.size()) };
    const int n_act { static_cast<int>(los.size() / dims::vec) };
    geo.lat.resize(n_alt * n_act);
    geo.lon.resize(n_alt * n_act);
    geo.height.resize(n_alt * n_act);
    geo.sza.resize(n_alt * n_act);
    geo.saa.resize(n_alt * n_act);
    geo.vza.resize(n_alt * n_act);
    geo.vaa.resize(n_alt * n_act);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < n_alt; ++i_alt) {
        // Compute the J2000-to-ECEF rotation quaternion and the sun
        // vector in ECEF.
        Quaternion q_j2000_ecef {};
        std::array<double, dims::vec> sun {};
        solarModel(tai_seconds[i_alt], tai_subsec[i_alt], sun, q_j2000_ecef);
        // Transform orbit positions from J2000 to ECEF
        std::array<double, dims::vec> orb_pos_ecef {
            orb_pos_j2000[i_alt * dims::vec + 0],
            orb_pos_j2000[i_alt * dims::vec + 1],
            orb_pos_j2000[i_alt * dims::vec + 2]
        };
        q_j2000_ecef.rotate(orb_pos_ecef);
        // Compute the spacecraft-to-ECEF attitude quaternions
        Quaternion q_sc_ecef { q_j2000_ecef };
        q_sc_ecef.multiply(att_quat_sc_j2000[i_alt]);
        for (int i_act {}; i_act < n_act; ++i_act) {
            std::array<double, dims::vec> los_cur {
                los[i_act * dims::vec + 0],
                los[i_act * dims::vec + 1],
                los[i_act * dims::vec + 2]
            };
            q_sc_ecef.rotate(los_cur);
            // LOS intersection point with ground
            std::array<double, dims::vec> pos {};
            intersectEllipsoid(los_cur, orb_pos_ecef, pos);
            // Current ALT and ACT index
            const int idx { i_alt * n_act + i_act };
            if (dem_filename.empty()) {
                cart2geo(pos, geo.lat[idx], geo.lon[idx], geo.height[idx]);
            } else {
                intersectTerrain(dem,
                                 los_cur,
                                 pos,
                                 geo.lat[idx],
                                 geo.lon[idx],
                                 geo.height[idx]);
            }
            solarAndViewingGeometry(sun,
                                    los_cur,
                                    geo.lat[idx],
                                    geo.lon[idx],
                                    geo.sza[idx],
                                    geo.saa[idx],
                                    geo.vza[idx],
                                    geo.vaa[idx]);
        }
    }
}

} // namespace tango
