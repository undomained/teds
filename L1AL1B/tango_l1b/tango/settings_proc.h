// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef SETTINGS_PROC_H
#define SETTINGS_PROC_H

#include "header.h"
#include "settings.h"

// Forward declaration.
class Settings_main;

// Collocations of switches for detector calibration steps.
struct Calibration_options {
    bool dark_apply { true };
    bool nonlin_apply { true };
    bool prnu_apply { true };
    bool rad_apply { true };

    // Optional calibration settings.
    uint32_t nonlin_niter = 100; // Number of iterations that will be allowed when performing the inversion that applies the non-linearity correction.
    double nonlin_tol = 1.0e-3; // Convergence tolerance when performing the inversion that applies the non-linearity correction.
    int stray_van_cittert_steps = 3; // Number of Van Cittert iterations to perform straylight correction.
};

// Calibration settings. This is a base class. Each CKD has her own settings
// class that inherit from this base class.
class Settings_proc : public Settings {

    // Protected constructor to prevent any raw Settings_proc instances to be created.
    protected:
    Settings_proc(
        Logger *creator
    );

    public:
    // Virtual destructor to prevent surprises.
    // If a subclass object is destructed when being a Settings pointer, the
    // destructor of the subclass is only called when the Settings object has a
    // virtual destructor. Note that the virtual destructor itself is also called
    // when a subclass has its own destructor.
    virtual ~Settings_proc() {}

    // Settings for everyone.
    bool skip = false; // Used for calibration steps that can be skipped.
    // (with ckd_file_in) and perform new steps. Only level 1+ output of current steps are saved.
    // Flag to skip viewports. Always legal, but only makes sense for viewport-related
    // CKD steps or L1B.
    vector<int> vp_skip; // Flags to skip (or destroy) viewports.
    vector<string> l1a_files; // Not used for logistic-only steps. For flight L1A, only one file.
    int l1a_type = L1A_CALIBRATION; // Flag for L1A manager style. Usage of this flag should be minimized.
    vector<uint64_t> image; // Single images to be processed. Only for flight L1A.
    vector<uint64_t> image_start; // Starts of image streaks to be processed. Only for flight L1A.
    vector<uint64_t> image_end; // Ends of image streaks to be processed. Only for flight L1A.

    // From main settings. You can adapt this per processor if you want.
    string binningtable_filename = ""; // File with all the used binning tables.
    vector<size_t> l1a_discard_frames_start; // Number of frames to discard at the beginning to mitigate memory effects.
    vector<size_t> l1a_discard_frames_end; // Number of frames to discard at the beginning to mitigate memory effects.
    // Omitting the l1a_discard_frames_start and l1a_discard_frames_end means
    // that you keep all the frames (e.g. for simulated measurements). If
    // your main settings sets the default to discard somw, but you want
    // to overwrite this to 'discard nothing' for one process, set both
    // to a single zero. This explicitly overrules the main settings.
    // Normally, omitting these settings in the processor triggers the
    // fallback to the defaults set in the main.

    // Interpret these arrays. This is done by the L1A manager.
    size_t l1a_discard_frames_attempts;

    int output_level = NC_FILL_INT; // Flag for additional detailed (non-operational) output.
    // Any optional output of level 1 or higher is discarded when reading previous CKD
    int output_level_max = 0; // Maximum supported output level. If set to zero, also no empty groups of detailed output are created. Set in constructor of derived class.

    // Flags for CKD usage.
    // These flags are there to use the same CKD slightly differently.
    // Examples are number of iterations for iterative calibration or a
    // convergence criterion. The idea is that these flags should stay
    // fixed operationally. To verify consistency, any data created will
    // have all these flags in their attributes.
    // TODO: Execute this write step. Where to put it into the L1B?
    // Maybe inside some structure.
    // It makes no sense to set these flags for CKD generation steps
    // that come before the step to which the flags are related.
    uint32_t &nonlin_niter = opt.nonlin_niter; // Number of iterations that will be allowed when performing the inversion that applies the non-linearity correction.
    double &nonlin_tol = opt.nonlin_tol; // Convergence tolerance when performing the inversion that applies the non-linearity correction.
    int &stray_van_cittert_steps = opt.stray_van_cittert_steps; // Number of Van Cittert iterations to perform straylight correction.
    bool& dark_apply { opt.dark_apply };
    bool& nonlin_apply { opt.nonlin_apply };
    bool& prnu_apply { opt.prnu_apply };
    bool& rad_apply { opt.rad_apply };

    // Apply main settings to the processor.
    int apply_main_settings(
        Settings_main *set_main // Main settings, which are already read.
    );

    // Interpreted contents. TODO: Maybe improve on this.
    bool overall_percentages = false; // Do a progress bar for the loop over batches. This is only nice for L1B as long there is no logging on screen inside the module.

    // L1A request flags, to be interpreted.
    bool l1a_navigation = false; // Request to read and interpolate navigation data (True only for L1B processing).
    // Possibly, we could add request flags for all the GSE-data.
    // Because the GSE-data is not final yet, we just read what we see.
    // Also, a decision has to be made for illumination levels for PRNU measurements.

    // L1X output.
    // There can be L1X output after any step. The user should provide
    // exactly the same number of L1X output files as there are L1A input
    // files and the L1A manager will write them in the same style as the
    // L1A file(s) (flight, calibration or frames).
    // Turn on a specific step of L1X output by defining their output files.
    // Define zero files for not choosing that output. Define exactly the
    // same number of files as L1A files to turn on that output.
    vector<string> l1x_outputfiles_raw; // Name of L1X output files after just co-adding and converting to floats.
    vector<string> l1x_outputfiles_dark; // Name of L1X output files after dark correction.
    vector<string> l1x_outputfiles_noise; // Name of L1X output files after noise characterization.
    vector<string> l1x_outputfiles_nonlin; // Name of L1X output files after non-linearity correction.
    vector<string> l1x_outputfiles_prnu; // Name of L1X output files after PRNU correction.
    vector<string> l1x_outputfiles_unbin; // Name of L1X output files after unbinning.
    vector<string> l1x_outputfiles_stray; // Name of L1X output files after straylight correction.
    vector<string> l1x_outputfiles_rebin; // Name of L1X output files after rebinning.
    vector<string> l1x_outputfiles_fov; // Name of L1X output files for extraced spectra in detector counts..
    vector<string> l1x_outputfiles_rad; // Name of L1X output files for extracted spectra in radiance units.

    // Detector calibration options.
    Calibration_options opt;

    protected:
    // Overwritten function for processor (not main) settings.
    int init_common(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) override;

    // Routine to overwrite, where step-specific settings are recognized.
    virtual int init_step(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) = 0;

};

#endif
