// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef NOISECAL_H
#define NOISECAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_noisecal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_noisecal(
        Logger *creator
    );
    ~Settings_noisecal();
    // Specific settings.
    size_t nfrac = 1; // Number of fractions to cut the L1A images.
    // Pixel mask criteria.
    double mask_g_min = NC_FILL_DOUBLE; // Minimum signal-dependent noise term.
    double mask_g_max = NC_FILL_DOUBLE; // Maximum signal-dependent noise term.
    double mask_n_min = NC_FILL_DOUBLE; // Minimum signal-independent noise term.
    double mask_n_max = NC_FILL_DOUBLE; // Maximum signal-independent noise term.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;
}; // }}}

class Noisecal : public Processor { // {{{

    // A public constructor.
    public:
    Noisecal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Noisecal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_noisecal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(
        size_t ibatch
    ) override;

    private:
    size_t nfile; // Number of L1A files (different illuminations).
    vector<size_t> nl1a_file; // Number of L1A (frames) per file (measurement).

    // Detailed output.
    // Dimension nfile is also used here.
    vector<double> det_signal; // Signal averages, dimensions (npix,nfile).
    vector<double> det_noise; // Signal standard deviaion, the square root of the variance, dimensions (npix,nfile).
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

}; // }}}

#endif
