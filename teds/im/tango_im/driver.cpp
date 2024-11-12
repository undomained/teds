// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver.h"

#include "forward_models.h"
#include "settings_im.h"

#include <netcdf>
#include <spdlog/spdlog.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/io.h>
#include <tango_l1b/l1.h>
#include <tango_l1b/timer.h>

namespace tango {

// Meta data such as the exposure time are read from the L1A input
// file but ignored by the instrument model. Instead, they are set by
// user settings. The L1A-L1B processor, however, only reads from the
// input file because after the IM run these parameter should be
// fixed.
static auto setL1Meta(const SettingsIM& settings,
                      std::vector<L1>& l1_products) -> void
{
    for (auto& l1 : l1_products) {
        l1.binning_table_id =
          static_cast<uint8_t>(settings.detector.binning_table_id);
        l1.nr_coadditions =
          static_cast<uint16_t>(settings.detector.nr_coadditions);
        l1.exposure_time = settings.detector.exposure_time;
    }
}

// Generate an ISRF kernel for convolving the radiances from the
// line-by-line grid onto the intermediate grids from the CKD swath
// section.
static auto genISRFKernel(const double fwhm_gauss,
                          const double shape,
                          const std::vector<double>& lbl_wavelengths,
                          const std::vector<double>& wavelengths,
                          Eigen::SparseMatrix<double>& isrf)
  -> void
{
    isrf =
      Eigen::SparseMatrix<double>(wavelengths.size(), lbl_wavelengths.size());
    const bool do_gaussian { std::abs(shape - 2.0) < 1e-30 };
    const double sigma { fwhm_gauss / (2.0 * std::sqrt(2.0 * std::log(2.0))) };
    const double sigma_2inv { 1.0 / (2.0 * sigma * sigma) };
    std::vector<Eigen::Triplet<double>> triplets {};
    std::vector<double> row_sums(isrf.rows(), 0.0); // For normalization
    for (int i_wave {}; i_wave < isrf.rows(); ++i_wave) {
        for (int i_lbl {}; i_lbl < isrf.cols(); ++i_lbl) {
            // Convolution result for this CKD wavelength
            double conv { wavelengths[i_wave] - lbl_wavelengths[i_lbl] };
            // Reduce the computational cost by considering the
            // limited extent of the Gaussian.
            constexpr double wave_rel_threshold { 1.5 };
            if (std::abs(conv) > wave_rel_threshold * fwhm_gauss) {
                continue;
            }
            if (do_gaussian) {
                conv = std::exp(-(conv * conv * sigma_2inv));
            } else {
                conv = std::pow(
                  2.0, -std::pow(2 * std::abs(conv) / fwhm_gauss, shape));
            }
            row_sums[i_wave] += conv;
            triplets.push_back({ i_wave, i_lbl, conv });
        }
    }
    isrf.setFromTriplets(triplets.begin(), triplets.end());
    // Normalize
    for (int k {}; k < isrf.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it { isrf, k }; it;
             ++it) {
            it.valueRef() /= row_sums[it.row()];
        }
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
    // For undoing nonlinearity calibration we need the inverse of the
    // nonlinearity spline from the CKD.
    const LinearSpline nonlin_spline { ckd.nonlin.spline.invert() };

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

    // Generate the ISRF matrix
    Eigen::SparseMatrix<double> isrf {};
    if (settings.isrf.enabled && l1_products.front().level == ProcLevel::sgm
        && settings.cal_level <= ProcLevel::sgm) {
        spdlog::info("Generating ISRF matrix for convolution");
        genISRFKernel(settings.isrf.fwhm_gauss,
                      settings.isrf.shape,
                      (*l1_products.front().wavelengths).front(),
                      ckd.swath.wavelengths,
                      isrf);
    }

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
            applyISRF(ckd, settings.isrf.enabled, isrf, l1);
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
            mapToDetector(ckd, settings.swath.b_spline_order, l1);
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
            nonlinearity(ckd, settings.nonlin.enabled, nonlin_spline, l1);
            timers[static_cast<int>(ProcLevel::nonlin)].stop();
        }
        // Dark current
        if (l1.level >= ProcLevel::dark_current
            && settings.cal_level < ProcLevel::dark_current) {
            timers[static_cast<int>(ProcLevel::dark_current)].start();
            darkCurrent(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark_current)].stop();
        }
        // Noise
        if (l1.level >= ProcLevel::noise
            && settings.cal_level < ProcLevel::noise) {
            timers[static_cast<int>(ProcLevel::noise)].start();
            noise(ckd, settings.noise.enabled, settings.noise.seed, l1);
            timers[static_cast<int>(ProcLevel::noise)].stop();
        }
        // Dark offset
        if (l1.level >= ProcLevel::dark_offset
            && settings.cal_level < ProcLevel::dark_offset) {
            timers[static_cast<int>(ProcLevel::dark_offset)].start();
            darkOffset(ckd, settings.dark.enabled, l1);
            timers[static_cast<int>(ProcLevel::dark_offset)].stop();
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
    if (l1_products.front().wavelengths->front().size()
        == ckd.wave.wavelengths.front().size()) {
        *l1_products.front().wavelengths = ckd.wave.wavelengths;
    } else {
        *l1_products.front().wavelengths =
          std::vector<std::vector<double>>(ckd.n_act, ckd.swath.wavelengths);
    }

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
    spdlog::info("         Dark current: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::dark_current)].time());
    spdlog::info("               Noise: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::noise)].time());
    spdlog::info("         Dark offset: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::dark_offset)].time());
    spdlog::info("      ADC conversion: {:8.3f} s",
                 timers[static_cast<int>(ProcLevel::l1a)].time());
    spdlog::info("      Writing output: {:8.3f} s", timer_output.time());
    spdlog::info("               Total: {:8.3f} s", timer_total.time());

    printHeading("Success");
}

} // namespace tango
