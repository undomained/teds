// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_main.h"
#include "settings_proc.h"

// Constructor that progresses the argument ultimately to its grandparent.
Settings_proc::Settings_proc(
    Logger *creator
) : Settings(creator)
{
}

// Read the settings file and stores the contents in the structure.
int Settings_proc::init_common( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{

    // Default settings, for every processor (or almost every processor).
    recognize_setting(skip); // Skip for calibration steps can can be skipped.
    recognize_setting_vector(vp_skip); // Flags to skip (or destroy) viewports.
    recognize_setting(output_level); // Flag for additional detailed (non-operational) output.
    recognize_setting_vector(l1a_files); // File or files of L1A input.
    recognize_setting(l1a_type); // L1A file type: Calibration or flight.
    recognize_setting_vector(image); // Single images to be processed.
    recognize_setting_vector(image_start); // Starts of image streaks to be processed.
    recognize_setting_vector(image_end); // Ends of image streaks to be processed.
    // Settings of which the default from the main settings can be overwritten.
    recognize_setting(binningtable_filename); // File with all the used binning tables.
    recognize_setting_vector(l1a_discard_frames_start); // Number of frames to discard at the beginning to mitigate memory effects.
    recognize_setting_vector(l1a_discard_frames_end); // Number of frames to discard at the beginning to mitigate memory effects.
    recognize_setting_vector(l1x_outputfiles_raw); // Name of L1X output files after just co-adding and converting to floats.
    recognize_setting_vector(l1x_outputfiles_dark); // Name of L1X output files after dark correction.
    recognize_setting_vector(l1x_outputfiles_noise); // Name of L1X output files after noise characterization.
    recognize_setting_vector(l1x_outputfiles_nonlin); // Name of L1X output files after non-linearity correction.
    recognize_setting_vector(l1x_outputfiles_prnu); // Name of L1X output files after PRNU correction.
    recognize_setting_vector(l1x_outputfiles_unbin); // Name of L1X output files after unbinning.
    recognize_setting_vector(l1x_outputfiles_stray); // Name of L1X output files after straylight correction.
    recognize_setting_vector(l1x_outputfiles_rebin); // Name of L1X output files after rebinning.
    recognize_setting_vector(l1x_outputfiles_fov); // Name of L1X output files for extracted spectra still in detector units.
    recognize_setting_vector(l1x_outputfiles_rad); // Name of L1X output files for extracted spectra in radiance units.

    // CKD usage settings.
    recognize_setting(nonlin_niter); // Number of iterations that will be allowed when performing the inversion that applies the non-linearity correction.
    recognize_setting(nonlin_tol); // Convergence tolerance when performing the inversion that applies the non-linearity correction.
    recognize_setting(stray_van_cittert_steps); // Number of Van Cittert iterations to perform straylight correction.
    recognize_setting(dark_apply);
    recognize_setting(nonlin_apply);
    recognize_setting(prnu_apply);
    recognize_setting(rad_apply);
    recognize_setting(stray_interpolating);

    // Own settings. This routine is overwritten by the inheriting class.
    handle(init_step(stream,key,value,recognized));

    return 0;

} // }}}

// Apply main settings that should be shared to all the processors.
// This is the case when a setting is actually relevant for processors,
// but the user wants to put it into the main section.
int Settings_proc::apply_main_settings( // {{{
    Settings_main *set_main
)
{
    // Copy relevant part of the main settings.
    if (binningtable_filename.compare("") == 0) binningtable_filename = set_main->binningtable_filename; // File with all the used binning tables.
    if (l1a_discard_frames_start.size() == 0) l1a_discard_frames_start = set_main->l1a_discard_frames_start; // Copy the vector.
    if (l1a_discard_frames_end.size() == 0) l1a_discard_frames_end = set_main->l1a_discard_frames_end; // Copy the vector.
    if (output_level == NC_FILL_INT) output_level = set_main->output_level; // Level of detailed output.

    return 0;

} // }}}

