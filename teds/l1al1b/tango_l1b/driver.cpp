// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver.h"

#include "binning_table.h"
#include "calibration.h"
#include "ckd.h"
#include "io.h"
#include "l1.h"
#include "settings_l1b.h"
#include "timer.h"

#include <spdlog/spdlog.h>

namespace tango {

auto driver(const SettingsL1B& settings,
            const int argc,
            const char* const argv[]) -> void
{
    // Set up loggers and print general information
    initLogging(false);
    printHeading("Tango L1B processor", false);
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
    CKD ckd(settings.io.ckd);

    // Initialize L1 products by reading all L1A data (everything is
    // stored in memory).
    std::vector<L1> l1_products {};
    readL1(settings.io.l1a,
           settings.image_start,
           settings.image_end.value_or(fill::i),
           l1_products);

    // Initialize the binning table and bin the CKD
    const BinningTable binning_table {
        ckd.n_detector_rows,
        ckd.n_detector_cols,
        settings.io.binning_table,
        static_cast<int>(l1_products.front().binning_table_id)
    };
    binning_table.bin(ckd.pixel_mask);
    binning_table.bin(ckd.dark.offset);
    binning_table.bin(ckd.dark.current);
    binning_table.bin(ckd.noise.g);
    binning_table.bin(ckd.noise.n2);
    binning_table.bin(ckd.prnu.prnu);

    // Run retrieval
    printHeading("Retrieval");
    std::array<Timer, static_cast<int>(ProcLevel::n_levels)> timers {};
    Timer timer_total {};
    timer_total.start();
#pragma omp parallel for schedule(dynamic)
    for (int i_alt = 0; i_alt < static_cast<int>(l1_products.size()); ++i_alt) {
        printPercentage(i_alt, l1_products.size(), "Processing images");
        auto& l1 { l1_products[i_alt] };
        // Normalize by bin sizes and unbin detector image if requested
        if (l1.level == ProcLevel::l1a
            && settings.cal_level >= ProcLevel::raw) {
            binningTable(ckd, binning_table, l1);
        }
        // Dark offset
        if (l1.level < ProcLevel::dark_offset
            && settings.cal_level >= ProcLevel::dark_offset) {
            timers[static_cast<int>(ProcLevel::dark_offset)].start();
            darkOffset(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark_offset)].stop();
        }
        // Noise
        if (l1.level < ProcLevel::noise
            && settings.cal_level >= ProcLevel::noise) {
            timers[static_cast<int>(ProcLevel::noise)].start();
            noise(ckd, settings.noise.enabled, binning_table, l1);
            timers[static_cast<int>(ProcLevel::noise)].stop();
        }
        // Dark current
        if (l1.level < ProcLevel::dark_current
            && settings.cal_level >= ProcLevel::dark_current) {
            timers[static_cast<int>(ProcLevel::dark_current)].start();
            darkCurrent(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark_current)].stop();
        }
        // Nonlinearity
        if (l1.level < ProcLevel::nonlin
            && settings.cal_level >= ProcLevel::nonlin) {
            timers[static_cast<int>(ProcLevel::nonlin)].start();
            nonlinearity(ckd, settings.nonlin.enabled, l1);
            timers[static_cast<int>(ProcLevel::nonlin)].stop();
        }
        // PRNU and QE
        if (l1.level < ProcLevel::prnu
            && settings.cal_level >= ProcLevel::prnu) {
            timers[static_cast<int>(ProcLevel::prnu)].start();
            prnu(ckd, settings.prnu.enabled, l1);
            timers[static_cast<int>(ProcLevel::prnu)].stop();
        }
        if (!l1.image.empty()) {
            removeBadValues(ckd, l1);
        }
        // Stray light
        if (l1.level < ProcLevel::stray
            && settings.cal_level >= ProcLevel::stray) {
            timers[static_cast<int>(ProcLevel::stray)].start();
            strayLight(
              ckd, binning_table, settings.stray.van_cittert_steps, l1);
            timers[static_cast<int>(ProcLevel::stray)].stop();
        }
        // Swath
        if (l1.level < ProcLevel::swath
            && settings.cal_level >= ProcLevel::swath) {
            timers[static_cast<int>(ProcLevel::swath)].start();
            mapFromDetector(
              ckd, binning_table, settings.swath.b_spline_order, l1);
            timers[static_cast<int>(ProcLevel::swath)].stop();
        }
        // Interpolate from intermediate to the main CKD wavelength
        // grids if necessary.
        if (!l1.spectra.empty()
            && l1.spectra.size() / ckd.n_act
                 != ckd.wave.wavelengths.front().size()) {
            changeWavelengthGrid(ckd, l1);
        }
        // Radiometric
        if (l1.level < ProcLevel::l1b && settings.cal_level >= ProcLevel::l1b) {
            timers[static_cast<int>(ProcLevel::l1b)].start();
            radiometric(ckd, settings.rad.enabled, l1);
            timers[static_cast<int>(ProcLevel::l1b)].stop();
        }
    }
    timer_total.stop();
    spdlog::info("Processing images 100.0%");
    // For writing to output, store the CKD wavelength grid in L1
    l1_products.front().wavelengths =
      std::make_shared<std::vector<std::vector<double>>>(ckd.wave.wavelengths);

    // Placeholder until we have geolocation
    copyGeometry(
      settings.io.l1a, settings.io.geometry, settings.image_start, l1_products);

    // Write output
    Timer timer_output {};
    timer_output.start();
    writeL1(settings.io.l1b, settings.getConfig(), l1_products, argc, argv);
    timer_output.stop();

    // Overview of timings
    spdlog::info("");
    spdlog::info("Timings:");
    spdlog::info("         Dark offset: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::dark_offset)].time());
    spdlog::info("               Noise: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::noise)].time());
    spdlog::info("        Dark current: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::dark_current)].time());
    spdlog::info("        Nonlinearity: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::nonlin)].time());
    spdlog::info("                PRNU: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::prnu)].time());
    spdlog::info("         Stray light: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::stray)].time());
    spdlog::info("               Swath: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::swath)].time());
    spdlog::info("         Radiometric: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::l1b)].time());
    spdlog::info("      Writing output: {:8.3f} s", timer_output.time());
    spdlog::info("               Total: {:8.3f} s", timer_total.time());

    printHeading("Success");
}

} // namespace tango
