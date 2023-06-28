// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef RADCAL_H
#define RADCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_radcal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_radcal(
        Logger *creator
    );
    ~Settings_radcal();
    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class Radcal : public Processor { // {{{

    // A public constructor.
    public:
    Radcal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Radcal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_radcal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(size_t ibatch, const Calibration_options& opt) override;

    private:
    // Batch-to-viewport coupling. Non-trivial if viewports are skipped.
    vector<size_t> ivp_batch;

}; // }}}

#endif
