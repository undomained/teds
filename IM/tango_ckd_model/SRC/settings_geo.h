#ifndef SETTINGS_GEO_H
#define SETTINGS_GEO_H

#include "header.h"
#include "logger.h"
#include "settings.h"

class Settings_geo : public Settings { // {{{

    public:
    // Constructor.
    Settings_geo(
        Logger *creator
    );
    // Specific settings.
    string utcfile = ""; // UTC Time-difference file.
    string demfile = ""; // Detailed elevation map. Leave empty for placeholder DEM function.
    double semi_major_axis = NC_FILL_DOUBLE; // Semi-major axis of the earth.
    double semi_minor_axis = NC_FILL_DOUBLE; // Semi-minor axis of the earth.
    double latitude_tol = NC_FILL_DOUBLE; // Tolerance for iterative calculation of latitude during geolocation (radians).
    double mountain = NC_FILL_DOUBLE; // Highest mountain to expect.
    double movedistance = NC_FILL_DOUBLE; // Distance to move in one iteration avoiding to skip mountains.
    double extremeweight = NC_FILL_DOUBLE; // Clip on weight factor for guess during geolocation convergence.
    double geolocation_tol = NC_FILL_DOUBLE; // Tolerance for convergence along line of sight (length units).

    // Overwritten virtual function.
    protected:
    int init_common(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

#endif
