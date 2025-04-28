// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver.h"

#include "forward_models.h"
#include "settings_im.h"

#include <common/ckd.h>
#include <common/io.h>
#include <common/isrf.h>
#include <common/l1.h>
#include <common/timer.h>
#include <spdlog/spdlog.h>

namespace tango {

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
                    settings.io_files.binning_table);

    // Read in the CKD
    printHeading("Reading CKD and input data");
    const CKD ckd { settings.io_files.ckd };
    // For undoing nonlinearity calibration we need the inverse of the
    // nonlinearity spline from the CKD.
    const LinearSpline nonlin_spline { ckd.nonlin.spline.invert() };

    // Read and initialize data
    L1 l1_prod {};
    readL1(settings.io_files.sgm,
           settings.alt_beg,
           settings.alt_end,
           l1_prod,
           settings.isrf.in_memory);
    l1_prod.exposure_time = settings.detector.exposure_time;

    // Read or construct the ISRF
    ISRF isrf {};
    if (settings.isrf.tabulated) {
        spdlog::info("Reading ISRF from file:");
        isrf.fromFile(
          settings.io_files.isrf, l1_prod.wavelengths, ckd.wave.wavelengths);
    } else {
        spdlog::info("Generating ISRF from generalized Gaussian parameters");
        isrf.fromGauss(l1_prod.wavelengths,
                       ckd.wave.wavelengths,
                       settings.isrf.fwhm,
                       settings.isrf.shape);
    }

    // Run the forward model
    printHeading("Forward model");
    Timer timer {};
    timer.start();
    // ISRF
    if (l1_prod.level == ProcLevel::sgm
        && settings.cal_level <= ProcLevel::sgm) {
        spdlog::info("ISRF convolution");
        applyISRF(
          ckd, settings.isrf.enabled, isrf, settings.io_files.sgm, l1_prod);
    }
    // Radiometric
    if (l1_prod.level >= ProcLevel::l1b
        && settings.cal_level < ProcLevel::l1b) {
        spdlog::info("Radiometric");
        radiometric(ckd, settings.rad.enabled, l1_prod);
    }
    // Swath (draw spectra on detector)
    if (l1_prod.level >= ProcLevel::swath
        && settings.cal_level < ProcLevel::swath) {
        spdlog::info("Detector mapping");
        mapToDetector(ckd, settings.swath.b_spline_order, l1_prod);
    }
    // Stray light
    if (l1_prod.level >= ProcLevel::stray
        && settings.cal_level < ProcLevel::stray) {
        spdlog::info("Stray light");
        strayLight(ckd, settings.stray.enabled, l1_prod);
    }
    // PRNU and QE
    if (l1_prod.level >= ProcLevel::prnu
        && settings.cal_level < ProcLevel::prnu) {
        spdlog::info("PRNU");
        prnu(ckd, settings.prnu.enabled, l1_prod);
    }
    // Nonlinearity
    if (l1_prod.level >= ProcLevel::nonlin
        && settings.cal_level < ProcLevel::nonlin) {
        spdlog::info("Nonlinearity");
        nonlinearity(settings.nonlin.enabled, nonlin_spline, l1_prod);
    }
    // Dark current
    if (l1_prod.level >= ProcLevel::dark_current
        && settings.cal_level < ProcLevel::dark_current) {
        spdlog::info("Dark signal");
        darkCurrent(ckd, settings.dark.enabled, l1_prod);
    }
    // Noise
    if (l1_prod.level >= ProcLevel::noise
        && settings.cal_level < ProcLevel::noise) {
        spdlog::info("Noise");
        // Only scale noise by coaddition factor if target level is
        // L1A. Otherwise the signal will not be scaled by coadditions
        // and neither should the noise be scaled.
        const int n_coadditions { settings.cal_level == ProcLevel::l1a
                                    ? settings.detector.nr_coadditions
                                    : 1 };
        noise(ckd,
              settings.noise.enabled,
              settings.noise.seed,
              n_coadditions,
              l1_prod);
    }
    // Dark offset
    if (l1_prod.level >= ProcLevel::dark_offset
        && settings.cal_level < ProcLevel::dark_offset) {
        spdlog::info("Dark offset");
        darkOffset(ckd, settings.dark.enabled, l1_prod);
    }
    // Optimal coadding factor and exposure time
    if (settings.optimal_coadd.enabled) {
        spdlog::info("Optimal coadding factors and exposure times:");
        estimateOptimalCoadd(ckd,
                             settings.optimal_coadd.FMC,
                             settings.detector.exposure_time,
                             settings.optimal_coadd.f_sat,
                             settings.optimal_coadd.full_well,
                             settings.optimal_coadd.t_dead,
                             l1_prod);
    }
    // Bin detector images if they exist
    if (l1_prod.signal.size() > 0) {
        spdlog::info("Detector image binning ({}x1)",
                     static_cast<int>(settings.detector.binning_table_id));
        binDetectorImages(ckd.n_detector_rows,
                          ckd.n_detector_cols,
                          settings.io_files.binning_table,
                          settings.detector.binning_table_id,
                          settings.cal_level > ProcLevel::l1a,
                          l1_prod);
    }
    // Coaddition and ADC conversion
    if (l1_prod.level >= ProcLevel::raw
        && settings.cal_level == ProcLevel::l1a) {
        spdlog::info("Analog-to-digital conversion");
        digitalToAnalog(settings.detector.nr_coadditions, l1_prod);
    }

    writeL1(settings.io_files.l1a,
            settings.getConfig(),
            l1_prod,
            settings.compress,
            argc,
            argv);

    if (!settings.io_files.navigation.empty()) {
        copyNavigationData(settings.io_files.navigation, settings.io_files.l1a);
    }

    timer.stop();
    spdlog::info("Total time: {:8.3f} s", timer.time());

    printHeading("Success");
}

} // namespace tango
