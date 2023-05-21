// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_main.h"

// Constructor: Set the tag to the right content.
Settings_main::Settings_main(
    Logger *creator
) : Settings(creator)
{
    tag = "main";
}
Settings_main::~Settings_main() {}
// Read the settings file and stores the contents in the structure.
int Settings_main::init_common( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{

    recognize_setting(ckd_file_in); // CKD input file (not for first step).
    recognize_setting(ckd_file_out); // CKD output file if there is one.
    recognize_setting_vector(process); // Processes to include.
    recognize_setting(last_calibration_step); // Process name of the last calibration step to execute (Processes need to support this, probably that will only be l1b).
    recognize_setting(binningtable_filename); // File with all the used binning tables.
    recognize_setting_vector(l1a_discard_frames_start); // Number of frames to discard at the beginning to mitigate memory effects.
    recognize_setting_vector(l1a_discard_frames_end); // Number of frames to discard at the beginning to mitigate memory effects.
    recognize_setting(output_level); // Flag for additional detailed (non-operational) output.

    return 0;
} // }}}

