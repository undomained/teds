#ifndef SETTINGS_MAIN_H
#define SETTINGS_MAIN_H

#include "header.h"
#include "settings.h"

// Main settings.
class Settings_main : public Settings {

    public:

    // Define settings here.
    int input_enhance = 1; // Factor to enhance the input spectral sampling to reduce numeric errors in ISRF convolution.
    bool photons = false; // Flag for input in photon units.
    string ckd_file_in = ""; // CKD input file (not for first step).
    string ckd_file_out = ""; // Not used.
    string binningtable_filename = ""; // File with all the used binning tables.
    int binning_table_id; // Binning table to use for the L1A.
    string l1x_input = ""; // L1X input file name.
    vector<uint64_t> image; // Single images to be processed.
    vector<uint64_t> image_start; // Starts of image streaks to be processed.
    vector<uint64_t> image_end; // Ends of image streaks to be processed.
    string l1a_outputfile = ""; // Name of L1A output file.
    string l1x_outputfile_raw = ""; // Name of L1X output files after just co-adding and converting to floats.
    string l1x_outputfile_dark = ""; // Name of L1X output files after dark correction.
    string l1x_outputfile_noise = ""; // Name of L1X output files after noise characterization.
    string l1x_outputfile_nonlin = ""; // Name of L1X output files after non-linearity correction.
    string l1x_outputfile_prnu = ""; // Name of L1X output files after PRNU correction.
    string l1x_outputfile_unbin = ""; // Name of L1X output files after ubninning correction (this is the same as PRNU, but in different structure).
    string l1x_outputfile_stray = ""; // Name of L1X output files after straylight correction.
    string l1x_outputfile_rebin = ""; // Name of L1X output files after rebinning (this is the same as STRAY, but in a different structure).
    string l1x_outputfile_fov = ""; // Name of L1X output files for extraced spectra in detector counts..
    string l1x_outputfile_rad = ""; // Name of L1X output files for extracted spectra in radiance units.
    double t_dwell { -1.0 };
    double full_well {};
    double f_sat {};
    double t_dead {};
    double exposure_time {};
    int nr_coadditions {};

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

};

#endif
