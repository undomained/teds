// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver.h"

#include "settings_im.h"
#include "uncalibration.h"

#include <netcdf>
#include <spdlog/spdlog.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/io.h>
#include <tango_l1b/l1.h>
#include <tango_l1b/timer.h>

namespace tango {

// Read the knots and values from the CKD file but construct the
// inverse of the nonlinearity spline.
static auto inverseNonlinearity(const std::string& ckd_file)
{
    const netCDF::NcFile nc { ckd_file, netCDF::NcFile::read };
    const auto grp { nc.getGroup("nonlinearity") };
    const auto n_knots { grp.getDim("knots").getSize() };
    std::vector<double> knots(n_knots);
    std::vector<double> y(n_knots);
    grp.getVar("knots").getVar(knots.data());
    grp.getVar("y").getVar(y.data());
    return LinearSpline { y, knots };
}

// Meta data such as the exposure time are read from the L1A input
// file but ignored by the instrument model. Instead, they are set by
// user settings. The L1A-L1B processor, however, only reads from the
// input file because after the IM run these parameter should be
// fixed.
static auto setL1Meta(const SettingsIM& settings, std::vector<L1>& l1_products)
  -> void
{
    for (auto& l1 : l1_products) {
        l1.binning_table_id =
          static_cast<uint8_t>(settings.detector.binning_table_id);
        l1.nr_coadditions =
          static_cast<uint16_t>(settings.detector.nr_coadditions);
        l1.exposure_time = settings.detector.exposure_time;
    }
}

auto driver(const SettingsIM& settings,
            const int argc,
            const char* const argv[]) -> void
{
    // Set up loggers and print general information
    initLogging(false);
    printHeading("Tango instrument model", false);
    printSystemInfo(TANGO_PROJECT_VERSION,
                    TANGO_GIT_COMMIT_ABBREV,
                    TANGO_CMAKE_HOST_SYSTEM,
                    TANGO_EXECUTABLE,
                    TANGO_CXX_COMPILER,
                    TANGO_CXX_COMPILER_FLAGS,
                    TANGO_LIBRARIES,
                    settings.io.binning_table);

    // Read in the CKD
    printHeading("Reading CKD and input data");
    const CKD ckd { settings.io.ckd };

    printHeading("Initialize the binning table");
    // Initialize the binning table
    const BinningTable binning_table { ckd.n_detector_rows,
                                       ckd.n_detector_cols,
                                       settings.io.binning_table,
                                       settings.detector.binning_table_id };

    // Read and initialize data
    std::vector<L1> l1_products {};
    readL1(settings.io.sgm,
           settings.image_start,
           settings.image_end.value_or(fill::i),
           l1_products);
    setL1Meta(settings, l1_products);

    // Run the forward model (main loop)
    printHeading("Forward model");
    std::array<Timer, static_cast<int>(ProcLevel::n_levels)> timers {};
    Timer timer_total {};
    timer_total.start();
#pragma omp parallel for schedule(dynamic)
    for (int i_alt = 0; i_alt < static_cast<int>(l1_products.size()); ++i_alt) {
        printPercentage(i_alt, l1_products.size(), "Processing images");
        auto& l1 { l1_products[i_alt] };

        // ISRF
        if (l1.level == ProcLevel::sgm
            && settings.cal_level <= ProcLevel::sgm) {
            timers[static_cast<int>(ProcLevel::sgm)].start();
            applyISRF(ckd, settings.isrf.enabled, settings.isrf.fwhm_gauss, l1);
            timers[static_cast<int>(ProcLevel::sgm)].stop();
        }
        // Radiometric
        if (l1.level >= ProcLevel::l1b && settings.cal_level < ProcLevel::l1b) {
            timers[static_cast<int>(ProcLevel::l1b)].start();
            radiometric(ckd, settings.rad.enabled, l1);
            timers[static_cast<int>(ProcLevel::l1b)].stop();
        }
        // Swath (draw on detector)
        if (l1.level >= ProcLevel::swath
            && settings.cal_level < ProcLevel::swath) {
            timers[static_cast<int>(ProcLevel::swath)].start();
            drawOnDetector(ckd, l1);
            timers[static_cast<int>(ProcLevel::swath)].stop();
        }
        // Stray light
        if (l1.level >= ProcLevel::stray
            && settings.cal_level < ProcLevel::stray) {
            timers[static_cast<int>(ProcLevel::stray)].start();
            strayLight(ckd, settings.stray.enabled, l1);
            timers[static_cast<int>(ProcLevel::stray)].stop();
        }
        // PRNU and QE
        if (l1.level >= ProcLevel::prnu
            && settings.cal_level < ProcLevel::prnu) {
            timers[static_cast<int>(ProcLevel::prnu)].start();
            prnu(ckd, settings.prnu.enabled, l1);
            timers[static_cast<int>(ProcLevel::prnu)].stop();
        }
        // Nonlinearity
        if (l1.level >= ProcLevel::nonlin
            && settings.cal_level < ProcLevel::nonlin) {
            timers[static_cast<int>(ProcLevel::nonlin)].start();
            // For undoing nonlinearity calibration we need the inverse of the
            // nonlinearity spline from the CKD.
            // Extra protection. Only do the inverseNonlinearity when non lin
            // ckd actually exists
            if (ckd.nonlin.enabled) {
                const LinearSpline nonlin_spline { inverseNonlinearity(
                  settings.io.ckd) };
                nonlinearity(ckd, settings.nonlin.enabled, nonlin_spline, l1);
            }
            timers[static_cast<int>(ProcLevel::nonlin)].stop();
        }
        // Dark current
        if (l1.level >= ProcLevel::dark
            && settings.cal_level < ProcLevel::dark) {
            timers[static_cast<int>(ProcLevel::dark)].start();
            darkCurrent(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark)].stop();
        }
        // Noise
        if (l1.level >= ProcLevel::noise
            && settings.cal_level < ProcLevel::noise) {
            timers[static_cast<int>(ProcLevel::noise)].start();
            noise(ckd, settings.noise.enabled, settings.noise.seed, l1);
            timers[static_cast<int>(ProcLevel::noise)].stop();
        }
        // Dark offset
        if (l1.level >= ProcLevel::dark
            && settings.cal_level < ProcLevel::dark) {
            timers[static_cast<int>(ProcLevel::dark)].start();
            darkOffset(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark)].stop();
        }
        // Coaddition, binning, and ADC conversion
        if (l1.level >= ProcLevel::raw
            && settings.cal_level == ProcLevel::l1a) {
            timers[static_cast<int>(ProcLevel::l1a)].start();
            digitalToAnalog(binning_table, l1);
            timers[static_cast<int>(ProcLevel::l1a)].stop();
        }
    }
    timer_total.stop();
    spdlog::info("Processing images 100.0%");
    // For writing to output switch out the LBL wavelength grid
    *l1_products.front().wavelength = ckd.wave.wavelength;

    Timer timer_output {};
    timer_output.start();
    writeL1(settings.io.l1a, settings.getConfig(), l1_products, argc, argv);
    timer_output.stop();

    spdlog::info("");
    spdlog::info("Timings:");
    spdlog::info("    ISRF convolution: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::sgm)].time());
    spdlog::info("         Radiometric: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::l1b)].time());
    spdlog::info("    Draw on detector: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::swath)].time());
    spdlog::info("         Stray light: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::stray)].time());
    spdlog::info("                PRNU: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::prnu)].time());
    spdlog::info("        Nonlinearity: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::nonlin)].time());
    spdlog::info("Dark offset, current: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::dark)].time());
    spdlog::info("               Noise: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::noise)].time());
    spdlog::info("      ADC conversion: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::l1a)].time());
    spdlog::info("      Writing output: {:8.3f} s", timer_output.time());
    spdlog::info("               Total: {:8.3f} s", timer_total.time());

    printHeading("Success");
}

} // namespace tango
