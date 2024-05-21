#pragma once

#include "header.h"
#include "settings.h"

namespace tango {

// ISRF settings.
class Settings_isrf : public Settings {

    public:

    // Define settings here.
    bool execute = false; // Execute ISRF convolution.
    double fwhm_gauss = NC_FILL_DOUBLE;

    // Constructor.
    Settings_isrf(
        Logger *creator
    );
    ~Settings_isrf(); // Destructor.

    // Overwritten common settings.
    protected:
    int init_common(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) override;

};

} // namespace tango
