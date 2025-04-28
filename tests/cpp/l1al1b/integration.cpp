// Integration tests for the L1A-L1B processor

#include "../testing.h"

#include <common/binning_table.h>
#include <common/io.h>
#include <l1al1b/calibration.h>
#include <l1al1b/driver.h>
#include <l1al1b/settings_l1b.h>

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
    settings.io_files.ckd = ckd_filename;
    settings.io_files.l1a = l1a_filename;
    settings.io_files.l1b = l1b_filename;
    settings.io_files.binning_table = binningtable_filename;
    settings.init();

    // Run all tests

    SECTION("Full chain, no binning")
    {
        // Run the simulator and read the L1B product from temporary space
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(l1.wavelengths.abs().sum(), WithinRel(163065.0, 1e-6));
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(1.6769428e21, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(6.9570764e18, 1e-6));
    }

    SECTION("Full chain, L1A binning 2")
    {
        // Bin L1A product and write to temporary space. The simulator
        // will read it from there and generate an L1B product.
        constexpr int bin_factor { 2 };
        const tango::BinningTable binning_table {
            0, 0, binningtable_filename, bin_factor
        };
        l1.signal = binning_table.binMulti(l1.signal, false);
        l1.binning_table_id = bin_factor;
        writeL1A(fixture_dir, l1a_filename, l1);
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(l1.wavelengths.abs().sum(), WithinRel(163065.0, 1e-6));
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(1.6933739e21, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(4.9733939e18, 1e-6));
    }

    SECTION("Full chain, L1B binning 5")
    {
        settings.bin_spectra = 5;
        tango::driver(settings);
        tango::readL1(l1b_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(l1.wavelengths.abs().sum(), WithinRel(163065.0, 1e-6));
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(3.3496472e20, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(6.2222817e17, 1e-6));
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
        CHECK_THAT(l1.wavelengths.abs().sum(), WithinRel(163065.0, 1e-6));
        CHECK_THAT(l1.spectra.abs().sum(), WithinRel(1.3641337e21, 1e-6));
        CHECK_THAT(l1.spectra_noise.abs().sum(), WithinRel(1.2243153e18, 1e-6));
        CHECK_THAT(geo.lat.abs().sum(), WithinRel(2590.0187939, 1e-6));
        CHECK_THAT(geo.lon.abs().sum(), WithinRel(722.7365457, 1e-6));
        CHECK_THAT(geo.height.abs().sum(), WithinAbs(2.42670388e-8, 1e-6));
        CHECK_THAT(geo.sza.abs().sum(), WithinRel(1959.7164232, 1e-6));
        CHECK_THAT(geo.saa.abs().sum(), WithinRel(7716.1682454, 1e-6));
        CHECK_THAT(geo.vza.abs().sum(), WithinRel(46.5414419, 1e-6));
        CHECK_THAT(geo.vaa.abs().sum(), WithinRel(4495.8437611, 1e-6));
    }

    // Teardown
    std::remove(ckd_filename.c_str());
    std::remove(l1a_filename.c_str());
    std::remove(l1b_filename.c_str());
    std::remove(binningtable_filename.c_str());
    std::remove(config_filename.c_str());
}
