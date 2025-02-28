// Integration tests for the L1A-L1B processor

#include "testing.h"

#include <tango_l1b/binning_table.h>
#include <tango_l1b/calibration.h>
#include <tango_l1b/driver.h>
#include <tango_l1b/io.h>
#include <tango_l1b/settings_l1b.h>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// Much of the interface is file-based
const std::string tmp_dir { std::filesystem::temp_directory_path() };
const std::string ckd_filename { tmp_dir + "/tango_ckd.nc" };
const std::string l1a_filename { tmp_dir + "/tango_l1a.nc" };
const std::string l1b_filename { tmp_dir + "/tango_l1b.nc" };
const std::string binningtable_filename { tmp_dir + "/tango_binningtable.nc" };
const std::string config_filename { tmp_dir + "/tango_config.yaml" };
const std::string fixture_dir { std::string(FIXTURE_DIR) };

TEST_CASE("integration tests")
{
    // Initialize data
    tango::CKD ckd {};
    readCKD(fixture_dir, ckd);
    writeCKD(fixture_dir, ckd_filename, ckd);

    tango::L1 l1 {};
    readL1A(fixture_dir, l1);
    writeL1A(fixture_dir, l1a_filename, l1);

    writeBinningTable(fixture_dir, binningtable_filename);

    // Initialize settings. For the settings class to work properly
    // the config should be read from a file even if it's mostly
    // empty.
    std::ofstream empty_config { config_filename };
    empty_config << "processing_version: test\n";
    empty_config.close();
    tango::SettingsL1B settings { config_filename };
    settings.io.ckd = ckd_filename;
    settings.io.l1a = l1a_filename;
    settings.io.l1b = l1b_filename;
    settings.io.binning_table = binningtable_filename;
    settings.init();

    // Run all tests

    SECTION("Full chain, no binning")
    {
        // Run the simulator and read the L1B product from temporary space
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(flatten2D(l1.wavelengths)),
                   WithinRel(5220516.0211424, 1e-6));
        CHECK_THAT(absSum(l1.spectra), WithinRel(1.0746193e21, 1e-6));
        CHECK_THAT(absSum(l1.spectra_noise), WithinRel(4.4664394e18, 1e-6));
    }

    SECTION("Full chain, L1A binning 2")
    {
        // Bin L1A product and write to temporary space. The simulator
        // will read it from there and generate an L1B product.
        constexpr int bin_factor { 2 };
        const tango::BinningTable binning_table {
            0, 0, binningtable_filename, bin_factor
        };
        binning_table.bin(l1.signal, false);
        l1.binning_table_id = bin_factor;
        writeL1A(fixture_dir, l1a_filename, l1);
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(flatten2D(l1.wavelengths)),
                   WithinRel(5220516.0211424, 1e-6));
        CHECK_THAT(absSum(l1.spectra), WithinRel(1.0852862e21, 1e-6));
        CHECK_THAT(absSum(l1.spectra_noise), WithinRel(3.1915396e18, 1e-6));
    }

    SECTION("Full chain, L1B binning 5")
    {
        settings.bin_spectra = 5;
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(flatten2D(l1.wavelengths)),
                   WithinRel(1044103.2042285, 1e-6));
        CHECK_THAT(absSum(l1.spectra), WithinRel(2.1492386e20, 1e-6));
        CHECK_THAT(absSum(l1.spectra_noise), WithinRel(8.9328787e17, 1e-6));
    }

    SECTION("Full chain, exact mapping")
    {
        settings.swath.exact_drawing = true;
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(flatten2D(l1.wavelengths)),
                   WithinRel(5220516.0211424, 1e-6));
        CHECK_THAT(absSum(l1.spectra), WithinRel(1.0791422e21, 1e-6));
        CHECK_THAT(absSum(l1.spectra_noise), WithinRel(4.1237624e18, 1e-6));
    }

    SECTION("Geolocation")
    {
        // Switch off most other processes
        settings.dark.enabled = false;
        settings.noise.enabled = false;
        settings.nonlin.enabled = false;
        settings.prnu.enabled = false;
        settings.stray.van_cittert_steps = 0;
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        tango::Geometry geo {};
        readGeo(l1b_filename, geo);
        CHECK_THAT(absSum(flatten2D(l1.wavelengths)),
                   WithinRel(5220516.0211424, 1e-6));
        CHECK_THAT(absSum(l1.spectra), WithinRel(8.7344179e20, 1e-6));
        CHECK_THAT(absSum(l1.spectra_noise), WithinRel(7.8356179e17, 1e-6));
        CHECK_THAT(absSum(geo.lat), WithinRel(2590.0187939, 1e-6));
        CHECK_THAT(absSum(geo.lon), WithinRel(722.7365457, 1e-6));
        CHECK_THAT(absSum(geo.height), WithinAbs(2.42670388e-8, 1e-6));
        CHECK_THAT(absSum(geo.sza), WithinRel(1959.7164232, 1e-6));
        CHECK_THAT(absSum(geo.saa), WithinRel(7716.1682454, 1e-6));
        CHECK_THAT(absSum(geo.vza), WithinRel(46.5414419, 1e-6));
        CHECK_THAT(absSum(geo.vaa), WithinRel(4495.8437611, 1e-6));
    }

    // Teardown
    std::remove(ckd_filename.c_str());
    std::remove(l1a_filename.c_str());
    std::remove(l1b_filename.c_str());
    std::remove(binningtable_filename.c_str());
    std::remove(config_filename.c_str());
}
