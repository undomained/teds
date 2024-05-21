// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef DIMCAL_H
#define DIMCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_dimcal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_dimcal(
        Logger *creator
    );
    ~Settings_dimcal();
    // Specific settings.
    size_t detector_spec = 0; // Number of detector pixels in the spectral dimension.
    size_t detector_spat = 0; // Number of detector pixels in the spatial dimension.
    size_t vp = 0; // Number of viewports.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;
}; // }}}

class Dimcal : public Processor { // {{{

    // A public constructor.
    public:
    Dimcal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Dimcal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_dimcal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;

}; // }}}

#endif
