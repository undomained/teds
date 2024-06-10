// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing all user defined configuration parameters of the
// Tango L1A-L1B processor. A single instance is used throughout the
// code.
//
// When adding a new parameter two places need to be updated: i)
// define the parameter by editing one of structs below or by defining
// a new struct and ii) call the scan function on the parameter in
// scanKeys so it can be read from the YAML file and/or printed on the
// screen. You may optionally add checks in Settings::checkParameters
// for things like correct file paths.

#pragma once

#include "settings.h"

namespace tango {

class SettingsL1B : public Settings
{
private:
    auto checkParameters() -> void override;

public:
    Setting<std::string> instrument {
        { "instrument" },
        "",
        "which instrument to run the simulator for (carbon, nitro)"
    };
    Setting<std::string> processing_version {
        { "processing_version" },
        {},
        "L0-L4 processing toolchain version"
    };
    Setting<ProcLevel> cal_level {
        { "cal_level" },
        ProcLevel::l1b,
        "If given, the last calibration step to be executed. Allowed values:\n"
        "dark, noise, nonlin, prnu, stray, swath, rad.",
    };
    Setting<bool> reverse_wavelength {
        { "reverse_wavelength" },
        false,
        "Whether to reverse the wavelength ordering in the L1B product. This\n"
        "affects wavelengths and observables such as radiance."
    };
    Setting<int> image_start {
        { "image_start" },
        {},
        "first image to be processed (counting starts at 0)",
    };
    Setting<std::optional<int>> image_end {
        { "image_end" },
        "last image to be processed (inclusive, counting starts at 0)",
    };

    struct
    {
        Setting<bool> enabled {
            { "dark", "enabled" },
            true,
            "whether to include dark offset and current calibration",
        };
    } dark;

    struct
    {
        Setting<bool> enabled {
            { "noise", "enabled" },
            true,
            "whether to include noise calibration",
        };
    } noise;

    struct
    {
        Setting<bool> enabled {
            { "nonlin", "enabled" },
            true,
            "whether to include nonlinearity calibration",
        };
    } nonlin;

    struct
    {
        Setting<bool> enabled {
            { "prnu", "enabled" },
            true,
            "whether to include PRNU calibration",
        };
    } prnu;

    struct
    {
        Setting<int> van_cittert_steps {
            { "stray", "van_cittert_steps" },
            3,
            "number of Van Cittert (deconvolution) steps for stray light\n"
            "correction (0 disables this step)",
        };
    } stray;

    struct
    {
        Setting<bool> enabled {
            { "swath", "enabled" },
            true,
            "whether to include geolocation",
        };
        Setting<double> spectrum_width {
            { "swath", "spectrum_width" },
            1.0,
            "width of each spectrum, in pixels, extracted from the detector"
        };
    } swath;

    struct
    {
        Setting<bool> enabled {
            { "rad", "enabled" },
            true,
            "whether to include radiometric calibration",
        };
    } rad;

    struct
    {
        Setting<std::string> ckd { { "io", "ckd" }, {}, "CKD file path" };
        Setting<std::string> binning_table {
            { "io", "binning_table" },
            {},
            "Path to a NetCDF file containing binning tables and count arrays\n"
            "Which ones to use is determined by the L1A fiel\n"
            "/image_attributes/binning_table"
        };
        Setting<std::string> l1a { { "io", "l1a" }, {}, "L1A product (input)" };
        Setting<std::string> l1b { { "io", "l1b" },
                                   {},
                                   "L1B product (output)" };
    } io;

    SettingsL1B() = default;
    SettingsL1B(const std::string& yaml_file) : Settings { yaml_file } {}
    auto scanKeys() -> void override;
    ~SettingsL1B() = default;
};

} // namespace tango
