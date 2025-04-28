// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Algorithms related to geometrical operations and geolocation

#pragma once

#include "constants.h"
#include "eigen.h"

namespace tango {

class DEM;

struct Geometry
{
    ArrayXXd lat {};
    ArrayXXd lon {};
    ArrayXXd height {};
    ArrayXXd vza {}; // Viewing zenith angle
    ArrayXXd vaa {}; // Viewing azimuth angle
    ArrayXXd sza {}; // Solar zenith angle
    ArrayXXd saa {}; // Solar azimuth angle
};

// Convert angles from radians to degrees. If dealing with longitudes,
// then first bring them to the -180/180 range.
auto moduloAndConvertAngles(ArrayXXd& angles) -> void;

// Given a series of timestamps and quaternions, find the quaternion
// corresponding to an arbitrary timestamp.
auto interpolateQuaternion(
  const Eigen::ArrayXd& att_times,
  const std::vector<Eigen::Quaterniond>& att_quaternions,
  const double image_time) -> Eigen::Quaterniond;

// Given a line-of-sight (LOS) vector and the satellite position in
// cartesian coordinates, return the point of intersection with the
// Earth ellipsoid in cartesian coordinates. Note that the LOS and
// satellite vector need to be rescaled which is why they are passed
// by value. Algorithm taken from Joint Polar Satellite System (JPPS)
// VIIRS Geolocation Algorithm Theoretical Basis Document (ATBD)
// section 3.3.2.2.1.
auto intersectEllipsoid(Eigen::Vector3d los,
                        Eigen::Vector3d sat) -> Eigen::Vector3d;

// Given a digital elevation model, a LOS vector, the satellite
// position in cartesian coordinates, and the ellipsoid intersection
// point computed by intersectEllipsoid, return the latitude,
// longitude, and height of the target point. Note that the
// intersection point is updated by this function.
auto intersectTerrain(DEM& dem,
                      const Eigen::Vector3d& los,
                      Eigen::Vector3d& pos,
                      double& lat,
                      double& lon,
                      double& height) -> void;

// Convert from cartesian to geodetic coordinates.
auto cart2geo(const Eigen::Vector3d& xyz,
              double& lat,
              double& lon,
              double& height) -> void;

// Given the sun and LOS vectors in ECR coordinates and the target
// point latitude and longitude, return the solar zenith and azimuth
// angles and the sensor viewing zenith and azimuth angles.
auto solarAndViewingGeometry(const Eigen::Vector3d& sun,
                             const Eigen::Vector3d& los,
                             const double lat,
                             const double lon,
                             double& sza,
                             double& saa,
                             double& vza,
                             double& vaa) -> void;

// Derive geolocation data from navigation data.
//  los - line-of-sight vectors of all ALT and ACT positions
//  tai_seconds - number of seconds since 1958-01-01
//  tai_subsec - fraction part of the number of seconds
//  orb_pos_j2000 - orbit positions in J2000
//  att_quat_sc_j2000 - spacecraft-to-J2000 attitude quaternions
//  geo - resulting viewing and solar geometries
auto geolocate(const std::string& dem_filename,
               const ArrayXNd<dims::vec>& los,
               const std::vector<uint32_t>& tai_seconds,
               const std::vector<double>& tai_subsec,
               const ArrayXNd<dims::vec>& orb_pos_j2000,
               const std::vector<Eigen::Quaterniond>& att_quat_sc_j2000,
               Geometry& geo) -> void;

} // namespace tango
