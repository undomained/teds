#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_isrf.h"

namespace tango {

// Constructor: Set the tag to the right content.
Settings_isrf::Settings_isrf(
    Logger *creator
) : Settings(creator)
{
    tag = "isrf";
}
Settings_isrf::~Settings_isrf() {}
// Read the settings file and stores the contents in the structure.
int Settings_isrf::init_common( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{
    recognize_setting(execute);
    recognize_setting(fwhm_gauss);

    return 0;

} // }}}

} // namespace tango
