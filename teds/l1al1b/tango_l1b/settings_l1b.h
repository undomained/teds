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
    Setting<std::string> processing_version {
        { "processing_version" },
        {},
        "L0-L4 processing toolchain version"
    };
    Setting<ProcLevel> cal_level {
        { "cal_level" },
        ProcLevel::l1b,
        "If given, the last calibration step to be executed. Allowed values:\n"
        "offset, noise, current, nonlin, prnu, stray, swath, l1b.",
    };
    Setting<size_t> alt_beg {
        { "alt_beg" },
        {},
        "first along-track position to be processed (counting starts at 0)",
    };
    Setting<std::optional<size_t>> alt_end {
        { "alt_end" },
        "last along-track position to be processed (inclusive, counting\n"
        "starts at 0)",
    };
    Setting<int> bin_spectra {
        { "bin_spectra" },
        1,
        "reduce the number of L1B spectra by averaging over this many spectra"
    };
    Setting<bool> compress { { "compress" },
                             true,
                             "whether to compress the L1B product" };

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
        Setting<int> b_spline_order {
            { "swath", "b_spline_order" },
            5,
            "order of 2D b-spline used for mapping spectra from the detector"
        };
        Setting<bool> geolocation {
            { "swath", "geolocation" },
            true,
            "Whether to include geolocation. If not then geometry is copied\n"
            "from the geometry file produced by the GM.",
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
        Setting<std::string> ckd { { "io_files", "ckd" }, {}, "CKD file path" };
        Setting<std::string> binning_table {
            { "io_files", "binning_table" },
            {},
            "Path to a NetCDF file containing binning tables and count arrays\n"
            "Which ones to use is determined by the L1A file\n"
            "/image_attributes/binning_table"
        };
        Setting<std::string> l1a {
            { "io_files", "l1a" },
            {},
            "L1A product (input). If cal_level, while running the instrument\n"
            "model, was set to anything other than l1a then this is actually\n"
            "not an L1A but a higher level product."
        };
        Setting<std::string> l1b {
            { "io_files", "l1b" },
            {},
            "L1B product (output). If cal_level is set to anything other than\n"
            "l1b then this is not an L1B but a lower level product."
        };
        Setting<std::string> dem { { "io_files", "dem" },
                                   {},
                                   "digital elevation model" };
        Setting<std::string> geometry { { "io_files", "geometry" },
                                        {},
                                        "geometry file (input)" };
    } io_files;

    SettingsL1B() = default;
    SettingsL1B(const std::string& yaml_file) : Settings { yaml_file } {}
    auto scanKeys() -> void override;
    ~SettingsL1B() = default;
};

} // namespace tango
