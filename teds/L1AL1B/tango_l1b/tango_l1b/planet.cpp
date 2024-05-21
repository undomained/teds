// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "vector.h"
#include "netcdf_object.h"
#include "planet.h"

namespace tango {

// Constructor.
Planet::Planet( // {{{
    Logger *creator,
    double semi_major_axis, // Semi-major axis [space unit].
    double semi_minor_axis, // Semi-minor aixs [space unit].
    double latitude_tolerance // Tolerance for searching latitude [radians].
) : Logger(creator)
{

    // Copy arguments to member variables.
    a = semi_major_axis;
    b = semi_minor_axis;
    tol = latitude_tolerance;

    // Square of eccentricity.
    ecc2 = 1.0 - pow(b,2.0) / pow(a,2.0);

} // }}}
Planet::~Planet(){}

void Planet::xyz( // {{{
    double lat, // Latitude.
    double lon, // Longitude.
    double alt, // Altitude.
    Vector *pos // Position (output).
)
{

    Vector &ref = *pos;
    ref[2] = a*(1.0-ecc2)*sin(lat) / sqrt(1.0-ecc2*pow(sin(lat),2.0)) + alt*sin(lat);
    double xy = a*cos(lat) / sqrt(1.0-ecc2*pow(sin(lat),2.0)) + alt*cos(lat);
    ref[0] = xy*cos(lon);
    ref[1] = xy*sin(lon);

} // }}}

void Planet::lla( // {{{
    Vector pos, // Position.
    double *lat, // Latitude (output).
    double *lon, // Longitude (output).
    double *alt // Altitude (output).
)
{

    // This is the inverse of xyz. It finds the latitude, longitude and altitude
    // from a point.
    // First find the point on the meridian between pole and pole with the correct
    // longitude. The longitude is never an issue.
    double xy = sqrt(pow(pos[0],2.0)+pow(pos[1],2.0)); // This one is positive.
    double z = pos[2];

    // Pole case (division by zero).
    if (xy == 0) {
        if (z>0) {
            *lat = 90.0*DEGREES;
            *alt = z-b;
        } else {
            *lat = -90.0*DEGREES;
            *alt = -z-b;
        }
        // Longitude not defined.
        *lon = 0.0;
        return;
    }

    // A lower limit of the latitude is the geocentric latitude.
    // Only if you choose an egg as earth semiminor axis greater than semimajor axis,
    // then, the geocentric latitude is an upper limit.
    // And the geodetic latitude of an inflated earth with the same eccentricity
    // is the other limit (I think).
    // An issue is that on the southern hemisphere, this is reversed, because they
    // are no upper and lower limits, but polar and equitorial limits.
    *lon = atan2(pos[1],pos[0]);
    double lim1 = atan(z/(xy*(1.0-ecc2)));
    double lim2;
    if (pow(xy,2.0)/pow(a,2.0) + pow(z,2.0)/pow(b,2.0) < 1.0) {
        double helpz = b * sqrt(1.0-pow(xy,2.0)/pow(a,2.0));
        if (z<0) helpz *= -1;
        lim2 = atan(helpz/(xy*(1.0-ecc2)));
    } else lim2 = atan(z/xy);
    vector<double> lats = {min(lim1,lim2),max(lim1,lim2)};
    while (lats[1]-lats[0] > tol) {
        double newlat = 0.5*(lats[0]+lats[1]);
        lats[northbound(newlat,xy,z)?0:1] = newlat;
    }
    *lat = 0.5*(lats[0]+lats[1]);
    // Calculate altitude adopting some stuff from the northbound function.
    double vec_xy = xy - a*cos(*lat) / sqrt(1.0-ecc2*pow(sin(*lat),2.0));
    double vec_z = z - a*(1.0-ecc2)*sin(*lat) / sqrt(1.0-ecc2*pow(sin(*lat),2.0));
    // But now, take the dot product with the zenith instead.
    double zen_xy = cos(*lat);
    double zen_z = sin(*lat);
    *alt = zen_xy*vec_xy + zen_z*vec_z;

} // }}}

bool Planet::northbound( // {{{
    double lat, // Attempted latitude.
    double xy, // Distance from pole axis.
    double z // Distance from equatorial plane.
)
{

    double vec_xy = xy - a*cos(lat) / sqrt(1.0-ecc2*pow(sin(lat),2.0));
    double vec_z = z - a*(1.0-ecc2)*sin(lat) / sqrt(1.0-ecc2*pow(sin(lat),2.0));
    double north_xy = -sin(lat);
    double north_z = cos(lat);
    // Dot product between zenith and north vector. This is a measure in which
    // the desire to go north is expressed.
    return vec_xy*north_xy + vec_z*north_z > 0.0;

} // }}}

int Planet::include_dem( // {{{
    string &filename
)
{

    // Open DEM file and read axes.
    nc_dem = make_unique<NetCDF_object>(this);
    handle(nc_dem->open(filename,NcFile::read));
    netcdf_check(nc_dem,dim_dem_lat = nc_dem->ncid->getDim("lat").getSize());
    netcdf_check(nc_dem,dim_dem_lon = nc_dem->ncid->getDim("lon").getSize());
    dem_lat.resize(dim_dem_lat);
    dem_lon.resize(dim_dem_lon);
    netcdf_check(nc_dem,nc_dem->ncid->getVar("lat").getVar(dem_lat.data()));
    netcdf_check(nc_dem,nc_dem->ncid->getVar("lon").getVar(dem_lon.data()));
    // Convert units.
    for (size_t ilat=0 ; ilat<dim_dem_lat ; ilat++) dem_lat[ilat] *= DEGREES;
    for (size_t ilon=0 ; ilon<dim_dem_lon ; ilon++) dem_lon[ilon] *= DEGREES;

    // The maps are not read because they are too large.
    netcdf_check(nc_dem,var_watermask = nc_dem->ncid->getVar("watermask"));
    netcdf_check(nc_dem,var_landheight = nc_dem->ncid->getVar("height"));
    netcdf_check(nc_dem,var_waterheight = nc_dem->ncid->getVar("water_surface_height"));

    return 0;

} // }}}

int Planet::dem( // {{{
    double lat, // Latitude in radians.
    double lon, // Longitude in radians.
    double &res // Resulting elevation.
)
{

    if (nc_dem) {
        // If there is a DEM NetCDF file, look up the DEM.
        // Search coordinates.

        // First forget about longitude wrapping.
        size_t ilon_left = NC_FILL_UINT64;
        size_t ilon_right = NC_FILL_UINT64;
        size_t ilat_left = NC_FILL_UINT64;
        size_t ilat_right = NC_FILL_UINT64;
        double weightleft_lat = NC_FILL_DOUBLE;
        double weightleft_lon = NC_FILL_DOUBLE;
        // If the latitude is out of range, we just take first/last value.
        if (lat < dem_lat[0]) {
            writelog(log_trace,"South pole extreme (test)");
            ilat_left = 0;
            ilat_right = 1;
            weightleft_lat = 1.0;
        } else if (lat > dem_lat[dim_dem_lat-1]) {
            writelog(log_trace,"North pole extreme (test)");
            ilat_left = dim_dem_lat - 2;
            ilat_right = dim_dem_lat - 1;
            weightleft_lat = 0.0;
        } else {
            ilat_left = 0;
            ilat_right = dim_dem_lat - 1;
            while (ilat_right-ilat_left != 1) {
                size_t ilat_new = (ilat_left+ilat_right)/2;
                if (lat > dem_lat[ilat_new]) ilat_left = ilat_new;
                else ilat_right = ilat_new;
            }
            weightleft_lat = (dem_lat[ilat_right]-lat) / (dem_lat[ilat_right] - dem_lat[ilat_left]);
        }
        // If the longitude is out of range, we will execute wrapping.
        bool wrap_lon = false;
        if (lon < dem_lon[0]) {
            writelog(log_trace,"West extreme (test)");
            ilon_left = dim_dem_lon - 1;
            ilon_right = 0;
            wrap_lon = true;
            weightleft_lon = (dem_lon[ilon_right] - lon) / (dem_lon[ilon_right] - dem_lon[ilon_left] + 360.0*DEGREES);
        } else if (lon > dem_lon[dim_dem_lon-1]) {
            writelog(log_trace,"East extreme (test)");
            ilon_left = dim_dem_lon - 1;
            ilon_right = 0;
            wrap_lon = true;
            weightleft_lon = (dem_lon[ilon_right] + 360.0*DEGREES - lon) / (dem_lon[ilon_right] - dem_lon[ilon_left] + 360.0*DEGREES);
        } else {
            ilon_left = 0;
            ilon_right = dim_dem_lon - 1;
            while (ilon_right-ilon_left != 1) {
                size_t ilon_new = (ilon_left+ilon_right)/2;
                if (lon > dem_lon[ilon_new]) ilon_left = ilon_new;
                else ilon_right = ilon_new;
            }
            weightleft_lon = (dem_lon[ilon_right]-lon) / (dem_lon[ilon_right] - dem_lon[ilon_left]);
        }
        // Verify that these variables are properly initialized.
        check_error(ilon_left == NC_FILL_UINT64 || ilon_right == NC_FILL_UINT64 || ilat_left == NC_FILL_UINT64 || ilat_right == NC_FILL_UINT64 || weightleft_lat == NC_FILL_DOUBLE || weightleft_lon == NC_FILL_DOUBLE,"Program error: Not all index or weights are properly initialized in index search algorithm.")
        vector<int8_t> watermask(4); // Char does not work, because NetCDF thinks it is a character. It is a one-byte integer.
        vector<int16_t> landheight(4);
        vector<int16_t> waterheight(4);
        if (wrap_lon) {
            // Read per longitude.
            netcdf_check(nc_dem,var_watermask.getVar({ilat_left,ilon_left},{2,1},watermask.data()));
            netcdf_check(nc_dem,var_watermask.getVar({ilat_left,ilon_right},{2,1},&watermask[2]));
            // Now, elements 1 and 2 should be switched.
            {
                char sw = watermask[1];
                watermask[2] = watermask[1];
                watermask[1] = sw;
            }
            netcdf_check(nc_dem,var_landheight.getVar({ilat_left,ilon_left},{2,1},landheight.data()));
            netcdf_check(nc_dem,var_landheight.getVar({ilat_left,ilon_right},{2,1},&landheight[2]));
            // Now, elements 1 and 2 should be switched.
            {
                int16_t sw = landheight[1];
                landheight[2] = landheight[1];
                landheight[1] = sw;
            }
            netcdf_check(nc_dem,var_waterheight.getVar({ilat_left,ilon_left},{2,1},waterheight.data()));
            netcdf_check(nc_dem,var_waterheight.getVar({ilat_left,ilon_right},{2,1},&waterheight[2]));
            // Now, elements 1 and 2 should be switched.
            {
                int16_t sw = waterheight[1];
                waterheight[2] = waterheight[1];
                waterheight[1] = sw;
            }
        } else {
            // Read in one go.
            netcdf_check(nc_dem,var_watermask.getVar({ilat_left,ilon_left},{2,2},watermask.data()));
            netcdf_check(nc_dem,var_landheight.getVar({ilat_left,ilon_left},{2,2},landheight.data()));
            netcdf_check(nc_dem,var_waterheight.getVar({ilat_left,ilon_left},{2,2},waterheight.data()));
        }
        vector<int16_t> dem_elev(4);
        for (size_t iel=0 ; iel<4 ; iel++) {
            if (watermask[iel] == 1) {
                // Water pixel. Use water surface height or zero if that
                // is unavailable.
                dem_elev[iel] = waterheight[iel];
                if (dem_elev[iel] == NC_FILL_SHORT) dem_elev[iel] = 0.0;
            } else {
                // Use land height. It cannot be fill value.
                dem_elev[iel] = landheight[iel];
            }
        }
        res = weightleft_lat*weightleft_lon*dem_elev[0] + weightleft_lat*(1.0-weightleft_lon)*dem_elev[1] + (1.0-weightleft_lat)*weightleft_lon*dem_elev[2] + (1.0-weightleft_lat)*(1.0-weightleft_lon)*dem_elev[3];
    } else {
        // No DEM NetCDF file, use placeholder DEM.
        res = 0.0;
    }

    return 0;

} // }}}

} // namespace tango
