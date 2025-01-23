// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing all user defined configuration parameters of the
// Tango L1A-L1B processor. A single instance is used throughout the
// code. For more details see settings.h and settings_l1b.h from the
// L1B processor.

#pragma once

#include <tango_l1b/settings.h>

namespace tango {

class SettingsIM : public Settings
{
private:
    auto checkParameters() -> void override;

public:
    Setting<std::string> instrument {
        { "instrument" },
        "",
        "which instrument to run the simulator for (Carbon, Nitro)"
    };
    Setting<std::string> processing_version {
        { "processing_version" },
        {},
        "L0-L4 processing toolchain version"
    };
    Setting<ProcLevel> cal_level {
        { "cal_level" },
        ProcLevel::l1a,
        "L1X level of the output file. For instance cal_level prnu means that\n"
        "the output data remains calibrated up to and including the PRNU\n"
        "correction, i.e. everything down to stray light has been unapplied."
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

    struct
    {
        Setting<int> binning_table_id {
            { "detector", "binning_table_id" },
            0,
            "which group from the binning table file to use",
        };
        Setting<double> exposure_time { { "detector", "exposure_time" },
                                        {},
                                        "exposure time for all images" };
        Setting<int> nr_coadditions {
            { "detector", "nr_coadditions" },
            1,
            "In the end, the images are multiplied by this number and the\n"
            "coadding factors are stored in the output file."
        };
    } detector;

    struct
    {
        Setting<bool> enabled {
            { "optimal_coadd", "enabled" },
            false,
            "whether to compute the optimal coadding factors and exposure times"
        };
        Setting<int> FMC {
            { "optimal_coadd", "FMC" },
            5,
            "Forward motion compensating factor for the optimal coadding\n"
            "estimation. Does not affect anything else about the simulation."
        };
        Setting<double> full_well { { "optimal_coadd", "full_well" },
                                    6679.665,
                                    "full well width" };
        Setting<double> t_dead { { "optimal_coadd", "t_dead" },
                                 0.006,
                                 "dead time, s" };
        Setting<double> f_sat { { "optimal_coadd", "f_sat" },
                                0.8,
                                "saturation threshold" };
    } optimal_coadd;

    struct
    {
        Setting<bool> enabled {
            { "isrf", "enabled" },
            true,
            "Whether to convolve the radiance spectra with the ISRF or\n"
            "linearly interpolate from the line-by-line onto the CKD\n"
            "wavelength grid.",
        };
        Setting<double> fwhm_gauss {
            { "isrf", "fwhm_gauss" },
            0.1,
            "the ISRF FWHM used for convolving the line-by-line spectra",
        };
        Setting<double> shape { { "isrf", "shape" },
                                2.0,
                                "ISRF shape parameter. Default is a Gauss." };
        Setting<bool> in_memory {
            { "isrf", "in_memory" },
            false,
            "Whether to load all spectra into memory before the convolution.\n"
            "Only applied when input is SGM spectra.",
        };
    } isrf;

    struct
    {
        Setting<bool> enabled {
            { "dark", "enabled" },
            true,
            "whether to include dark offset and current",
        };
    } dark;

    struct
    {
        Setting<bool> enabled {
            { "noise", "enabled" },
            true,
            "whether to include noise",
        };
        Setting<int> seed {
            { "noise", "seed" },
            0,
            "random number generator seed if applying noise",
        };
        // If swath.exact_drawing is true then do not bin the detector
        // image but instead artifically scale noise. Binning table ID
        // is always set to 1 in that case. The value is determined by
        // binning and the detector mapping algorithm choice and is
        // not a user parameter.
        double artificial_scaling { 1.0 };
    } noise;

    struct
    {
        Setting<bool> enabled {
            { "nonlin", "enabled" },
            true,
            "whether to include nonlinearity",
        };
    } nonlin;

    struct
    {
        Setting<bool> enabled {
            { "prnu", "enabled" },
            true,
            "whether to include PRNU",
        };
    } prnu;

    struct
    {
        Setting<bool> enabled {
            { "stray", "enabled" },
            true,
            "whether to include stray light",
        };

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
            "order of 2D b-spline used for mapping spectra to the detector"
        };
        Setting<bool> exact_drawing {
            { "swath", "exact_drawing" },
            false,
            "Whether to use the 'exact' drawing algorithm. In principle, it\n"
            "doesn't have any interpolation error but it yields an\n"
            "unphysical detector image and does not work with a large number\n"
            "of ACT angles."
        };
    } swath;

    struct
    {
        Setting<bool> enabled {
            { "rad", "enabled" },
            true,
            "whether to include radiometric",
        };
    } rad;

    struct
    {
        Setting<std::string> ckd { { "io", "ckd" }, {}, "CKD file path" };
        Setting<std::string> binning_table {
            { "io", "binning_table" },
            {},
            "Path to a NetCDF file containing binning tables and count arrays\n"
            "Which ones to use is determined by the L1A file variable\n"
            "/image_attributes/binning_table"
        };
        Setting<std::string> sgm {
            { "io", "sgm" },
            {},
            "Radiation scene from the scene generation module (SGM input).\n"
            "One can also use an L1B or a lower level product. The instrument\n"
            "then skips some of the first processin steps."
        };
        Setting<std::string> l1a {
            { "io", "l1a" },
            {},
            "L1A product (output). If cal_level was set to anything other\n"
            "than l1a then this is actually not an L1A but a higher level\n"
            "product."
        };
        Setting<std::string> geometry { { "io", "geometry" },
                                        {},
                                        "geometry file (input)" };
        Setting<std::string> navigation { { "io", "navigation" },
                                          {},
                                          "navigation data (input)" };
    } io;

    SettingsIM() = default;
    SettingsIM(const std::string& yaml_file) : Settings { yaml_file } {}
    auto scanKeys() -> void override;
    ~SettingsIM() = default;
};

} // namespace tango
