#include "header.h"
#include "settings.h"
#include "settings_geo.h"

namespace tango {

// Settings functions.
Settings_geo::Settings_geo( // {{{
    Logger *creator
) : Settings(creator)
{
    tag = "geo";
} // }}}

int Settings_geo::init_common( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(utcfile); // UTC Time-difference file.
    recognize_setting(demfile); // Detailed elevation map. Leave empty for placeholder DEM function.
    recognize_setting(semi_major_axis); // Semi-major axis of the earth.
    recognize_setting(semi_minor_axis); // Semi-minor axis of the earth.
    recognize_setting(latitude_tol); // Tolerance for iterative calculation of latitude during geolocation (radians).
    recognize_setting(mountain); // Highest mountain to expect.
    recognize_setting(movedistance); // Distance to move in one iteration avoiding to skip mountains.
    recognize_setting(extremeweight); // Clip on weight factor for guess during geolocation convergence.
    recognize_setting(geolocation_tol); // Tolerance for convergence along line of sight (length units).

    return 0;

} // }}}

} // namespace tango
