// Unit tests for the L1A-L1B processor

#include "testing.h"

#include <tango_l1b/binning_table.h>
#include <tango_l1b/calibration.h>

using Catch::Matchers::WithinRel;

TEST_CASE("unit tests")
{
    // Initialize data

    tango::CKD ckd {};
    readCKD(std::string(FIXTURE_DIR), ckd);

    tango::L1 l1 {};
    readL1(std::string(FIXTURE_DIR), l1);

    const tango::BinningTable binning_table {
        ckd.n_detector_rows, ckd.n_detector_cols, "", 0
    };

    // Run all tests

    SECTION("Dark offset")
    {
        tango::darkOffset(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(7627694.4368265271, 1e-6));
    }

    SECTION("Noise")
    {
        tango::noise(ckd, true, binning_table, 1.0, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(24981.7265996802, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8238421.598286096, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(10834706.9656610023, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(4563.8899392562, 1e-6));
    }

    SECTION("Stray light, 1 kernel")
    {
        tango::strayLight(ckd, binning_table, 3, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(9178760.4724600203, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }

    SECTION("Stray light, 2 kernels")
    {
        ckd.stray.n_kernels = 1;
        ckd.stray.kernel_rows.resize(ckd.stray.n_kernels);
        ckd.stray.kernel_cols.resize(ckd.stray.n_kernels);
        ckd.stray.kernels_fft.resize(ckd.stray.n_kernels);
        ckd.stray.kernel_fft_sizes.resize(ckd.stray.n_kernels);
        ckd.stray.weights.resize(ckd.stray.n_kernels);
        ckd.stray.edges.resize(ckd.stray.n_kernels * tango::box::n);
        tango::strayLight(ckd, binning_table, 3, l1);
        CHECK_THAT(absSum(l1.signal), WithinRel(9178937.454984488, 1e-6));
        CHECK_THAT(absSum(l1.noise), WithinRel(3264.0, 1e-6));
    }
}
