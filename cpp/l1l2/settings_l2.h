// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing all configuration parameters of the L2 processor

#pragma once

#include <common/settings.h>

namespace tango {

class SettingsL2 : public Settings
{
private:
    auto checkParameters() -> void override;

public:
    Setting<std::string> processing_version {
        { "processing_version" },
        {},
        "L0-L4 processing toolchain version"
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
    Setting<bool> compress { { "compress" },
                             true,
                             "whether to compress the L2 product" };

    struct
    {
        Setting<int> max_iter { { "retrieval", "max_iter" },
                                35,
                                "Gauss-Newton iteration limit" };
        Setting<double> chi2_lim {
            { "retrieval", "chi2_lim" },
            0.05,
            "chi2 convergence criterion between two Gauss-Newton iterations"
        };
        Setting<int> n_albedos { { "retrieval", "n_albedos" },
                                 2,
                                 "number of albedo coefficients to fit for" };
        Setting<std::vector<std::string>> gases { { "retrieval", "gases" },
                                                  { "CO2", "CH4", "H2O" },
                                                  "gases to be retrieved" };
        Setting<std::vector<double>> initial_concentrations {
            { "retrieval", "initial_concentrations" },
            { 400e-6, 1700e-9, 9000e-6 },
            "Initial concentrations of gases to be retrieved. Order must\n"
            "match those in [retrieval][gases]."
        };
    } retrieval;

    struct
    {
        Setting<int> n_layers { { "atmosphere", "n_layers" }, 30, "" };
        Setting<double> layer_thickness { { "atmosphere", "layer_thickness" },
                                          1000.0,
                                          "" };
        Setting<double> surface_pressure { { "atmosphere", "surface_pressure" },
                                           101300.0,
                                           "" };
    } atmosphere;

    struct
    {
        Setting<double> wave_start {
            { "spec_settings", "wave_start" },
            1590.0,
            "lower bound of L2 wavelength window for fitting, nm"
        };
        Setting<double> wave_end {
            { "spec_settings", "wave_end" },
            1675.0,
            "upper bound of L2 wavelength window for fitting, nm"
        };
        Setting<double> wave_extend {
            { "spec_settings", "wave_extend" },
            2.0,
            "L2 wavelength window extension for convolutions, nm"
        };
        Setting<double> dwave { { "spec_settings", "dwave" },
                                0.002,
                                "L2 wavelength window step size, nm" };
    } spec_settings;

    struct
    {
        Setting<bool> tabulated {
            { "isrf", "tabulated" },
            true,
            "whether to use the tabulated ISRF from [io_files][isrf] or\n"
            "generate it from the generalized Gaussian parameters"
        };
        Setting<double> fwhm {
            { "isrf", "fwhm" },
            0.456,
            "the ISRF FWHM used for convolving the line-by-line spectra",
        };
        Setting<double> shape { { "isrf", "shape" },
                                2.91,
                                "ISRF shape parameter. Default is a Gauss." };
    } isrf;

    struct
    {
        Setting<std::string> isrf {
            { "io_files", "isrf" },
            {},
            "ISRF, only needed [isrf][tabulated] is true, which is the default"
        };
        Setting<std::string> l1b { { "io_files", "l1b" },
                                   {},
                                   "L1B data product (input)" };
        Setting<std::string> atmosphere { { "io_files", "atmosphere" },
                                          {},
                                          "standard atmosphere file" };
        Setting<std::string> afgl { { "io_files", "afgl" }, {}, "" };
        Setting<std::string> sun_reference { { "io_files", "sun_reference" },
                                             {},
                                             "solar irradiance file" };
        Setting<std::string> dump_xsec {
            { "io_files", "dump_xsec" },
            {},
            "absorption cross-sections of the trace gases"
        };
        Setting<std::string> l2 { { "io_files", "l2" },
                                  {},
                                  "L2 data product (output)" };
        Setting<std::string> l2_diag {
            { "io_files", "l2_diag" },
            {},
            "File with extra L2 diagnostics information, i.e. those not\n"
            "already included in the main product (output)."
        };
    } io_files;

    SettingsL2() = default;
    SettingsL2(const std::string& yaml_file) : Settings { yaml_file } {}
    auto scanKeys() -> void override;
    ~SettingsL2() = default;
};

} // namespace tango
