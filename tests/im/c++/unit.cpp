// Unit tests for the instrument model

#include "../../l1al1b/c++/testing.h"

#include <tango_im/forward_models.h>

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
        CHECK_THAT(absSum(l1.image), WithinRel(8894711.5631734822, 1e-6));
    }

    SECTION("Noise")
    {
        tango::noise(ckd, true, 0, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(8261612.4975394309, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(3264.0, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(8283984.4017139031, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(ckd, true, ckd.nonlin.spline.invert(), l1);
        CHECK_THAT(absSum(l1.image), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(absSum(l1.image), WithinRel(6529884.2487379303, 1e-6));
        CHECK_THAT(absSum(l1.stdev), WithinRel(3264.0, 1e-6));
    }
}
