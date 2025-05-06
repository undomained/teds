// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver_l2.h"

#include "atmosphere.h"
#include "gauss_newton.h"
#include "optic_abs_prop.h"
#include "read_sun_spectrum.h"
#include "ref_profiles.h"
#include "settings_l2.h"

#include <common/ckd.h>
#include <common/io.h>
#include <common/isrf.h>
#include <common/l1.h>
#include <common/l2.h>
#include <common/timer.h>
#include <netcdf>
#include <numbers>
#include <spdlog/spdlog.h>

namespace tango {

auto driverL2(const SettingsL2& settings,
              const int argc,
              const char* const argv[]) -> void
{
    // Set up loggers and print general information
    initLogging();
    printHeading("Tango L2 processor", false);
    printSystemInfo(TANGO_PROJECT_VERSION,
                    TANGO_GIT_COMMIT_ABBREV,
                    TANGO_CMAKE_HOST_SYSTEM,
                    TANGO_EXECUTABLE,
                    TANGO_CXX_COMPILER,
                    TANGO_CXX_COMPILER_FLAGS,
                    TANGO_LIBRARIES);

    // Read in data
    printHeading("Reading input data");
    L1 l1b {};
    readL1(settings.io_files.l1b, settings.alt_beg, settings.alt_end, l1b);

    // Define line-by-line spectral grid
    const double& w_start { settings.spec_settings.wave_start };
    const double& w_end { settings.spec_settings.wave_end };
    const double& w_extend { settings.spec_settings.wave_extend };
    const double& dwave { settings.spec_settings.dwave };
    const int n_lbl { static_cast<int>((w_end - w_start + 2 * w_extend)
                                       / dwave) };
    const Eigen::ArrayXd wave_lbl =
      w_start - w_extend
      + dwave * Eigen::ArrayXd::LinSpaced(n_lbl, 0, n_lbl - 1);

    // Trim L1B data to match the L2 wavelength window
    const auto i_start { std::distance(
      l1b.wavelengths.begin(),
      std::ranges::find_if(
        l1b.wavelengths, [w_start](const double x) { return x > w_start; })) };
    const auto i_end { std::distance(
      l1b.wavelengths.begin(),
      std::ranges::find_if(l1b.wavelengths,
                           [w_end](const double x) { return x > w_end; })) };
    l1b.wavelengths = l1b.wavelengths(Eigen::seq(i_start, i_end)).eval();
    l1b.spectra = l1b.spectra(Eigen::all, Eigen::seq(i_start, i_end)).eval();
    l1b.spectra_noise =
      l1b.spectra_noise(Eigen::all, Eigen::seq(i_start, i_end)).eval();

    // Read atmosphere
    const Atmosphere atm { settings.atmosphere.n_layers,
                           settings.atmosphere.layer_thickness,
                           settings.atmosphere.surface_pressure,
                           settings.io_files.afgl };

    // Read optical absorption properties
    OpticAbsProp optics { settings.io_files.dump_xsec };
    // OpticAbsProp optics { settings.io_files.dump_xsec, wave_lbl };
    optics.checkWavelengthGrid(wave_lbl);
    optics.setOptDepthRef(atm, settings.retrieval.gases);

    // Define isrf function
    ISRF isrf {};
    if (settings.isrf.tabulated) {
        spdlog::info("Reading ISRF from file:");
        isrf.fromFile(settings.io_files.isrf, wave_lbl, l1b.wavelengths);
    } else {
        spdlog::info("Generating ISRF from generalized Gaussian parameters");
        isrf.fromGauss(
          wave_lbl, l1b.wavelengths, settings.isrf.fwhm, settings.isrf.shape);
    }

    // Solar irradiance: line-by-line and convolved
    Eigen::ArrayXd sun_lbl =
      readSunSpectrum(settings.io_files.sun_reference, wave_lbl);
    Eigen::ArrayXd sun = isrf.convolve(sun_lbl);

    // Need reference column profiles for fitting
    RefProfiles ref_profiles {};
    for (int i {}; i < static_cast<int>(settings.retrieval.gases.size()); ++i) {
        const std::string& gas { settings.retrieval.gases[i] };
        ref_profiles.gases[gas] = atm.gases.at(gas);
        ref_profiles.gas_sums[gas] = atm.gases.at(gas).sum();
        ref_profiles.initial[gas] =
          settings.retrieval.initial_concentrations[i];
    }

    // Initialize L2 product
    L2 l2 { l1b.n_alt,
            static_cast<int>(l1b.spectra.rows() / l1b.n_alt),
            static_cast<int>(l1b.wavelengths.size()),
            settings.atmosphere.n_layers,
            settings.retrieval.gases };

    // Timings

    // Run the forward model
    printHeading("Retrieval");
    Timer timer {};
    timer.start();

#pragma omp parallel for schedule(dynamic) firstprivate(optics)
    for (int i_gp = 0; i_gp < l1b.spectra.rows(); ++i_gp) {
        printPercentage(i_gp, l1b.spectra.rows(), "Processing ground points:");
        const double& sza { l1b.geo.sza.reshaped<Eigen::RowMajor>()(i_gp) };
        const double& vza { l1b.geo.vza.reshaped<Eigen::RowMajor>()(i_gp) };
        // Variances
        const Eigen::ArrayXd S_meas =
          l1b.spectra_noise.row(i_gp) * l1b.spectra_noise.row(i_gp);

        gaussNewton(i_gp,
                    settings,
                    ref_profiles,
                    atm,
                    wave_lbl,
                    sun_lbl,
                    l1b.spectra.row(i_gp),
                    S_meas,
                    std::cos(sza),
                    std::cos(vza),
                    isrf,
                    optics,
                    sun,
                    l2);

        if (!l2.converged(i_gp)) {
            // Decompose the ground point index into ALT and ACT indices
            const auto i_alt { i_gp / l2.converged.cols() };
            const auto i_act { i_gp - i_alt * l2.converged.cols() };
            spdlog::warn("Pixel ({},{}) not converged", i_alt, i_act);
        }
    }
    spdlog::info("Processing ground points: 100.00%");

    // Define proxy products
    constexpr double xch4_model { 1.8e-6 };
    constexpr double xco2_model { 410e-6 };
    l2.proxys.at("CO2") =
      l2.mixing_ratios.at("CO2") / l2.mixing_ratios.at("CH4") * xch4_model;
    l2.proxys.at("CH4") =
      l2.mixing_ratios.at("CH4") / l2.mixing_ratios.at("CO2") * xco2_model;
    const Eigen::ArrayXd rel_error =
      ((l2.precisions.at("CO2") / l2.mixing_ratios.at("CO2")).square()
       / (l2.precisions.at("CH4") / l2.mixing_ratios.at("CH4")).square())
        .sqrt();
    l2.proxy_precisions.at("CO2") = rel_error * l2.mixing_ratios.at("CO2");
    l2.proxy_precisions.at("CH4") = rel_error * l2.mixing_ratios.at("CH4");

    spdlog::info("");
    spdlog::info("Writing output");
    writeL2(settings.io_files.l2,
            settings.getConfig(),
            l2,
            l1b.geo,
            atm.z_lay,
            ref_profiles.gases,
            settings.compress,
            argc,
            argv);

    if (!settings.io_files.l2_diag.empty()) {
        writeL2Diagnostics(settings.io_files.l2_diag,
                           settings.getConfig(),
                           l1b,
                           l2,
                           settings.compress,
                           argc,
                           argv);
    }

    timer.stop();
    spdlog::info("");
    spdlog::info("Total time: {:8.3f} s", timer.time());

    printHeading("Success");
}

} // namespace tango
