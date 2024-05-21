#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_main.h"

namespace tango {

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

    // Put recognize_setting here.
    recognize_setting(input_enhance); // Factor to enhance the input spectral sampling to reduce numeric errors in ISRF convolution.
    recognize_setting(photons); // Flag for input in photon units.
    recognize_setting(log_file_path); // log file path .
    recognize_setting(ckd_file_in); // CKD input file (not for first step).
    recognize_setting(binningtable_filename); // File with all the used binning tables.
    recognize_setting(binning_table_id); // Binning table to use for the L1A.
    recognize_setting(l1x_input); // L1X input file name.
    recognize_setting_vector(image); // Single images to be processed.
    recognize_setting_vector(image_start); // Starts of image streaks to be processed.
    recognize_setting_vector(image_end); // Ends of image streaks to be processed.
    recognize_setting(l1a_outputfile); // Name of L1A output file.
    recognize_setting(l1x_outputfile_raw); // Name of L1X output files after just co-adding and converting to floats.
    recognize_setting(l1x_outputfile_dark); // Name of L1X output files after dark correction.
    recognize_setting(l1x_outputfile_noise); // Name of L1X output files after noise characterization.
    recognize_setting(l1x_outputfile_nonlin); // Name of L1X output files after non-linearity correction.
    recognize_setting(l1x_outputfile_prnu); // Name of L1X output files after PRNU correction.
    recognize_setting(l1x_outputfile_unbin); // Name of L1X output files after ubninning correction (this is the same as PRNU, but in different structure).
    recognize_setting(l1x_outputfile_stray); // Name of L1X output files after straylight correction.
    recognize_setting(l1x_outputfile_rebin); // Name of L1X output files after rebinning (this is the same as STRAY, but in a different structure).
    recognize_setting(l1x_outputfile_fov); // Name of L1X output files for extraced spectra in detector counts..
    recognize_setting(l1x_outputfile_rad); // Name of L1X output files for extracted spectra in radiance units.
    recognize_setting(t_dwell);
    recognize_setting(full_well);
    recognize_setting(f_sat);
    recognize_setting(t_dead);
    recognize_setting(exposure_time);
    recognize_setting(nr_coadditions);
    recognize_setting(dark_apply);
    recognize_setting(nonlin_apply);
    recognize_setting(stray_apply);
    recognize_setting(prnu_apply);
    recognize_setting(rad_apply);
    recognize_setting(swath_apply);
    recognize_setting(noise_apply);

    return 0;
} // }}}

} // namespace tango
