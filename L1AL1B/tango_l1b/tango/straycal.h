// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef STRAYCAL_H
#define STRAYCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

class CKD;

class Settings_straycal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_straycal (
        Logger *creator
    );
    ~Settings_straycal();

    // Specific settings.
    // When convolving with a moving kernel, the image is transformed
    // into a different size. The dimension of the new array is
    // 2*padding pixels larger than the original. We use -1 here as a
    // dummy default value. The actual default is
    // 0.1*ckd->dim_detector_spat.
    int padding { -1 };
    // Intensity threshold for pixels that make up a moving kernel. It
    // is a ratio of maximum intensity whre the total kernel has been
    // subtracted from a calibrations measurement. Anything less is
    // considered noise.
    double moving_kernel_intensity_threshold { 0.001 };
    bool dry_run { false };

    // Overwritten virtual function.
    protected:
    int init_step (
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class Straycal : public Processor { // {{{

    // A public constructor.
    public:
    Straycal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Straycal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_straycal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_batch(size_t ibatch, const Calibration_options& opt) override;
}; // }}}

#endif
