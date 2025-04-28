// Unit tests for the L1A-L1B processor

#include "../testing.h"

#include <common/binning_table.h>
#include <l1al1b/calibration.h>

using Catch::Matchers::WithinRel;

TEST_CASE("unit tests")
{
    // Initialize data

    tango::CKD ckd {};
    readCKD(std::string(FIXTURE_DIR), ckd);

    tango::L1 l1 {};
    readL1A(std::string(FIXTURE_DIR), l1);

    const tango::BinningTable binning_table {
        ckd.n_detector_rows, ckd.n_detector_cols, "", 0
    };

    // Run all tests

    SECTION("Dark offset")
    {
        tango::darkOffset(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(7600086.3545328, 1e-6));
    }

    SECTION("Noise")
    {
        tango::noise(ckd, true, binning_table, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(26884.6936471, 1e-6));
    }

    SECTION("Dark current")
    {
        tango::darkCurrent(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(26064709.4484070, 1e-6));
    }

    SECTION("Nonlinearity")
    {
        tango::nonlinearity(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8261203.0, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("PRNU")
    {
        tango::prnu(ckd, true, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(11443847.1976984, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(4624.8039625, 1e-6));
    }

    SECTION("Stray light, 2 kernels")
    {
        tango::strayLight(ckd, binning_table, 3, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(9178760.4724600203, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Stray light, 1 kernel")
    {
        ckd.stray.n_kernels = 1;
        ckd.stray.kernel_rows.conservativeResize(ckd.stray.n_kernels);
        ckd.stray.kernel_cols.conservativeResize(ckd.stray.n_kernels);
        ckd.stray.kernels_fft.resize(ckd.stray.n_kernels);
        ckd.stray.kernel_fft_sizes.conservativeResize(ckd.stray.n_kernels);
        ckd.stray.weights.conservativeResize(ckd.stray.n_kernels,
                                             ckd.stray.weights.cols());
        ckd.stray.edges.conservativeResize(ckd.stray.n_kernels, tango::box::n);
        tango::removeBadValues(ckd, l1);
        tango::strayLight(ckd, binning_table, 3, l1);
        CHECK_THAT(l1.signal.abs().sum(), WithinRel(8058139.1881834455, 1e-6));
        CHECK_THAT(l1.noise.abs().sum(), WithinRel(3264.0, 1e-6));
    }

    SECTION("Detector mapping")
    {
        tango::removeBadValues(ckd, l1);
        tango::mapFromDetector(ckd, binning_table, 3, l1);
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(11177936.6464352, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(5000.0, 1e-6));
    }

    SECTION("Radiometric")
    {
        tango::removeBadValues(ckd, l1);
        tango::mapFromDetector(ckd, binning_table, 3, l1);
        tango::radiometric(ckd, true, l1);
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(1.0241783e21, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(4.5812494e17, 1e-6));
    }
}
