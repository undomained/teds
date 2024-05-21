#pragma once

#include "header.h"
#include "settings.h"

namespace tango {

// Noise settings.
class Settings_noise : public Settings {

    public:

    bool noise_apply { false };
    unsigned long seed {};

    // Constructor.
    Settings_noise(
        Logger *creator
    );
    ~Settings_noise(); // Destructor.

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
