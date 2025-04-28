// Unit tests for the instrument model

#include "../testing.h"

#include <im/forward_models.h>

using Catch::Matchers::WithinRel;

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
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8922319.6454672, 1e-6));
    }

    SECTION("Noise")
    {
        tango::noise(ckd, true, 0, 1, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8261612.4975394309, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(26437548.8992291, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(true, ckd.nonlin.spline.invert(), l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(6175154.6307808431, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Stray light")
    {
        tango::strayLight(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(7435369.4285430796, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Detector mapping")
    {
        tango::mapToDetector(ckd, 3, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(2.7706936e20, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Radiometric")
    {
        tango::radiometric(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }
}
