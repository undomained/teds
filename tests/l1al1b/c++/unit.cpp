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

    // Run all tests

    SECTION("Dark offset")
    {
        tango::darkOffset(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(7600086.3545328248, 1e-6));
    }

    SECTION("Noise")
    {
        const tango::BinningTable binning_table {
            ckd.n_detector_rows, ckd.n_detector_cols, "", 0
        };
        tango::noise(ckd, true, binning_table, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(26884.6936471068, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(26064709.4484070353, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(11443847.1976984218, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(4624.8039624599, 1e-6));
    }
}
