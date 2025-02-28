// Unit tests for the instrument model

#include "../../l1al1b/c++/testing.h"

#include <tango_im/forward_models.h>

using Catch::Matchers::WithinRel;

// Use cubic splines to extend the spectra to match the intermediate
// wavelength grid dimension.
auto extendSpectra(const tango::CKD& ckd, tango::L1& l1) -> void
{
    const int n_waves_ext { static_cast<int>(ckd.swath.wavelengths.size()) };
    std::vector<double> spectra_ext(ckd.n_act * n_waves_ext);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        std::vector<double> x_values(ckd.n_detector_cols);
        std::vector<double> y_values(ckd.n_detector_cols);
        for (int i {}; i < ckd.n_detector_cols; ++i) {
            x_values[ckd.n_detector_cols - i - 1] =
              ckd.wave.wavelengths[i_act][i];
            y_values[ckd.n_detector_cols - i - 1] =
              l1.spectra[i_act * ckd.n_detector_cols + i];
        }
        const tango::CubicSpline spline { x_values, y_values };
        for (int j {}; j < n_waves_ext; ++j) {
            spectra_ext[i_act * n_waves_ext + j] =
              spline.eval(ckd.swath.wavelengths[j]);
        }
        l1.wavelengths[i_act] = ckd.swath.wavelengths;
    }
    l1.spectra = std::move(spectra_ext);
}

TEST_CASE("unit tests")
{
    // Initialize data

    tango::CKD ckd {};
    readCKD(std::string(FIXTURE_DIR), ckd);

    tango::L1 l1 {};
    readL1A(std::string(FIXTURE_DIR), l1);

    // Run all tests

    SECTION("Dark offset")
    {
        tango::darkOffset(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8894711.5631734822, 1e-6));
    }

    SECTION("Noise")
    {
        tango::noise(ckd, true, 0, 1, 1.0, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8261612.4975394309, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8283984.4017139031, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(ckd, true, ckd.nonlin.spline.invert(), l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(6529884.2487379303, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("Stray light")
    {
        tango::strayLight(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(7435369.4285430796, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("Detector mapping")
    {
        extendSpectra(ckd, l1);
        l1.spectra_noise.resize(l1.spectra.size());
        tango::mapToDetector(ckd, 3, false, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(2.9450128267645482e20, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("Radiometric")
    {
        extendSpectra(ckd, l1);
        l1.spectra_noise.resize(l1.spectra.size());
        tango::radiometric(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }
}
