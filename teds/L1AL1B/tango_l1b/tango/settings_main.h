// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef SETTINGS_MAIN_H
#define SETTINGS_MAIN_H

#include "header.h"
#include "settings.h"

// Main framework settings. These are there to launch different processors in
// a batch.
class Settings_main : public Settings {

    public:
    string ckd_file_in; // Not used for very first (DIM) step.
    string ckd_file_out; // Always used, except for L1B.
    vector<string> process; // Chosen processes to execute in string form.
    string last_calibration_step; // Process name of the last calibration step to execute (Processes need to support this, probably that will only be l1b).
    // Sets defaults to the processors. Then, the settings only have to be
    // given once. Processors can always overwrite if they want to.
    string binningtable_filename = ""; // File with all the used binning tables.
    string log_file_path = ""; // log file path.
    vector<size_t> l1a_discard_frames_start; // Number of frames to discard at the beginning to mitigate memory effects.
    vector<size_t> l1a_discard_frames_end; // Number of frames to discard at the beginning to mitigate memory effects.
    int output_level = 0; // Flag for additional detailed (non-operational) output.

    // Constructor.
    Settings_main(
        Logger *creator
    );
    ~Settings_main(); // Destructor.

    // Overwritten common settings.
    protected:
    int init_common(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) override;

    // No step-specific settings, because Settings_main has no children.

};

#endif
