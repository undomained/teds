#include "header.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_l1b.h"
#include "ckd.h"

namespace tango {

// Constructor: Set the tag to the right content.
Settings_l1b::Settings_l1b(
    Logger *creator
) : Settings(creator)
{
    tag = "l1b";
}
Settings_l1b::~Settings_l1b() {}
// Read the settings file and stores the contents in the structure.
int Settings_l1b::init_common( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{

    // Put recognize_setting here.
    recognize_setting(execute); // Flag for writing truth L1B file.
    recognize_setting(outputfile); // Output L1B file.
    recognize_setting_vector(l1b_wavelength); // Central wavelengths of polarimetric product.
    recognize_setting_vector(intensity_wavelength); // Wavelengths for intensity-only grid.
    recognize_setting(resolving_power); // Full-width half-maximum of Gauss convolution of polarimetric L1B product.
    recognize_setting(gauss_range); // Range of Gauss in FWHMs on both sides.

    return 0;

} // }}}

} // namespace tango
