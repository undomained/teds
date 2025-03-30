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
                    settings.io_files.binning_table);

    // Read in the CKD
    printHeading("Reading CKD and input data");
    CKD ckd(settings.io_files.ckd);

    // Initialize L1 products by reading all L1A data (everything is
    // stored in memory).
    L1 l1_prod {};
    readL1(settings.io_files.l1a, settings.alt_beg, settings.alt_end, l1_prod);

    // Initialize the binning table and bin the CKD
    const BinningTable binning_table { ckd.n_detector_rows,
                                       ckd.n_detector_cols,
                                       settings.io_files.binning_table,
                                       static_cast<int>(
                                         l1_prod.binning_table_id) };
    ckd.bin(binning_table);

    // Run retrieval
    printHeading("Retrieval");
    Timer timer {};
    timer.start();
    // Normalize detector image by bin sizes
    if (l1_prod.level == ProcLevel::l1a
        && settings.cal_level > ProcLevel::l1a) {
        spdlog::info("Scaling with bin size and coaddition factors");
        binScaling(ckd, binning_table, l1_prod);
    }
    // Dark offset
    if (l1_prod.level < ProcLevel::dark_offset
        && settings.cal_level >= ProcLevel::dark_offset) {
        spdlog::info("Dark offset");
        darkOffset(ckd, settings.dark.enabled, l1_prod);
    }
    // Noise
    if (l1_prod.level < ProcLevel::noise
        && settings.cal_level >= ProcLevel::noise) {
        spdlog::info("Noise");
        noise(ckd, settings.noise.enabled, binning_table, l1_prod);
    }
    // Dark current
    if (l1_prod.level < ProcLevel::dark_current
        && settings.cal_level >= ProcLevel::dark_current) {
        spdlog::info("Dark signal");
        darkCurrent(ckd, settings.dark.enabled, l1_prod);
    }
    // Nonlinearity
    if (l1_prod.level < ProcLevel::nonlin
        && settings.cal_level >= ProcLevel::nonlin) {
        spdlog::info("Nonlinearity");
        nonlinearity(ckd, settings.nonlin.enabled, l1_prod);
    }
    // PRNU and QE
    if (l1_prod.level < ProcLevel::prnu
        && settings.cal_level >= ProcLevel::prnu) {
        spdlog::info("PRNU");
        prnu(ckd, settings.prnu.enabled, l1_prod);
    }
    if (!l1_prod.signal.empty() && settings.cal_level >= ProcLevel::stray) {
        spdlog::info("Smoothing out bad detector signals");
        removeBadValues(ckd, l1_prod);
    }
    // Stray light
    if (l1_prod.level < ProcLevel::stray
        && settings.cal_level >= ProcLevel::stray) {
        spdlog::info("Stray light");
        strayLight(
          ckd, binning_table, settings.stray.van_cittert_steps, l1_prod);
    }
    // Swath
    if (l1_prod.level < ProcLevel::swath
        && settings.cal_level >= ProcLevel::swath) {
        spdlog::info("Detector mapping");
        mapFromDetector(
          ckd, binning_table, settings.swath.b_spline_order, l1_prod);
    }
    // Radiometric
    if (l1_prod.level < ProcLevel::l1b
        && settings.cal_level >= ProcLevel::l1b) {
        spdlog::info("Radiometric");
        radiometric(ckd, settings.rad.enabled, l1_prod);
    }
    if (settings.swath.geolocation) {
        spdlog::info("Geolocation");
        geolocate(settings.io_files.dem,
                  ckd.swath.los,
                  l1_prod.tai_seconds,
                  l1_prod.tai_subsec,
                  l1_prod.orb_pos,
                  l1_prod.att_quat,
                  l1_prod.geo);
    } else {
        spdlog::info("Copying geometry from geometry file");
        copyGeometry(settings.io_files.l1a,
                     settings.io_files.geometry,
                     settings.alt_beg,
                     l1_prod);
    }
    // Bin L1B spectra and geometry
    if (l1_prod.level >= ProcLevel::swath) {
        binL1B(settings.bin_spectra, l1_prod);
    }

    writeL1(settings.io_files.l1b, settings.getConfig(), l1_prod, argc, argv);

    timer.stop();
    spdlog::info("Total time: {:8.3f} s", timer.time());

    printHeading("Success");
}

} // namespace tango
