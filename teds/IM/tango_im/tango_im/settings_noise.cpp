#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_noise.h"

namespace tango {

// Constructor: Set the tag to the right content.
Settings_noise::Settings_noise(
    Logger *creator
) : Settings(creator)
{
    tag = "noise";
}
Settings_noise::~Settings_noise() {}
// Read the settings file and stores the contents in the structure.
int Settings_noise::init_common( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{
    recognize_setting(noise_apply);
    recognize_setting(seed);
    return 0;
} // }}}

 } // namespace tango
