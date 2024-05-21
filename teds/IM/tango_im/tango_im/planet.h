// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

// Forward declaration.
class Vector;
class NetCDF_object;

class Planet : public Logger {

    private:
    double a; // Semi-major axis.
    double b; // Semi-minor axis.
    double tol; // Tolerance for searching latitude [radians].
    double ecc2; // Eccentricity squared.
    // For DEM.
    unique_ptr<NetCDF_object> nc_dem;
    size_t dim_dem_lat; // Number of latitudes of detailed elevation map.
    size_t dim_dem_lon; // Number of longitudes of detailed elevation map.
    vector<double> dem_lat; // Latitude axis of detailed elevation map.
    vector<double> dem_lon; // Longitude axis of detailed elevation map.
    // An issue is that the elevation map shows the surface height or
    // the sea bottom (negative) height. For sea pixels, we need the
    // water surface height. The water surface height is fillvalue for
    // the land pixels, so a mask is needed.
    NcVar var_watermask; // Variable identifier for land/sea mask.
    NcVar var_landheight; // Variable identifier for land surface height.
    NcVar var_waterheight; // Variable identifier for water surface height.

    public:
    // Constructor.
    Planet(
        Logger *creator,
        double semi_major_axis, // Semi-major axis [space unit].
        double semi_minor_axis, // Semi-minor aixs [space unit].
        double latitude_tolerance // Tolerance for searching latitude [radians].
    );
    ~Planet();

    void xyz(
        double lat, // Latitude.
        double lon, // Longitude.
        double alt, // Altitude.
        Vector *pos // Position (output).
    );

    void lla(
        Vector pos, // Position.
        double *lat, // Latitude (output).
        double *lon, // Longitude (output).
        double *alt // Altitude (output).
    );

    int include_dem(
        string &filename // Filename for detailed elevation map.
    );

    // Acquisition of DEM, from file or using placeholder.
    int dem(
        double lat, // Latitude in radians.
        double lon, // Longitude in radians.
        double &res // Resulting elevation.
    );

    private:
    bool northbound(
        double lat, // Attempted latitude.
        double xy, // Distance from pole axis.
        double z // Distance from equatorial plane.
    );

};

} // namespace tango
