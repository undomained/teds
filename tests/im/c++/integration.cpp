// Integration tests for the instrument model

#include "../../l1al1b/c++/testing.h"

#include <tango_im/driver.h>
#include <tango_im/settings_im.h>
#include <tango_l1b/io.h>

#include <cstdio>

using Catch::Matchers::WithinRel;

// Much of the interface is file-based
const std::string tmp_dir { std::filesystem::temp_directory_path() };
const std::string ckd_filename { tmp_dir + "/tango_ckd.nc" };
const std::string sgm_filename { tmp_dir + "/tango_sgm.nc" };
const std::string l1a_filename { tmp_dir + "/tango_l1a.nc" };
const std::string geometry_filename { tmp_dir + "/tango_geometry.nc" };
const std::string navigation_filename { tmp_dir + "/tango_navigation.nc" };
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
    writeSGM(fixture_dir, sgm_filename, ckd, l1);

    writeGeometry(geometry_filename, l1);

    writeNavigation(fixture_dir, navigation_filename);

    writeBinningTable(fixture_dir, binningtable_filename);

    // Initialize settings. For the settings class to work properly
    // the config should be read from a file even if it's mostly
    // empty.
    std::ofstream empty_config { config_filename };
    empty_config << "processing_version: test\n";
    empty_config.close();
    tango::SettingsIM settings { config_filename };
    settings.detector.exposure_time = 0.01724385;
    settings.detector.nr_coadditions = 2;
    settings.isrf.tabulated = false;
    settings.isrf.in_memory = true;
    settings.isrf.fwhm = 0.5;
    settings.io_files.ckd = ckd_filename;
    settings.io_files.sgm = sgm_filename;
    settings.io_files.l1a = l1a_filename;
    settings.io_files.geometry = geometry_filename;
    settings.io_files.navigation = navigation_filename;
    settings.io_files.binning_table = binningtable_filename;
    settings.init();

    // Run all tests

    SECTION("Full chain")
    {
        // Run the simulator and read the L1A product from temporary space
        tango::driver(settings);
        tango::readL1(l1a_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(l1.signal), WithinRel(1284800.0, 1e-6));
    }

    SECTION("Full chain, no ADC or binning")
    {
        settings.cal_level = tango::ProcLevel::raw;
        tango::driver(settings);
        tango::readL1(l1a_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(l1.signal), WithinRel(642164.3174249, 1e-6));
    }

    SECTION("Full chain, binning 4")
    {
        settings.detector.binning_table_id = 4;
        tango::driver(settings);
        tango::readL1(l1a_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(l1.signal), WithinRel(1284464.0, 1e-6));
    }

    SECTION("Full chain, binning 4, no ADC")
    {
        settings.cal_level = tango::ProcLevel::raw;
        settings.detector.binning_table_id = 4;
        tango::driver(settings);
        tango::readL1(l1a_filename, 0, std::optional<size_t> {}, l1, true);
        CHECK_THAT(absSum(l1.signal), WithinRel(163950.7084508, 1e-6));
    }

    // Teardown
    std::remove(ckd_filename.c_str());
    std::remove(sgm_filename.c_str());
    std::remove(l1a_filename.c_str());
    std::remove(geometry_filename.c_str());
    std::remove(navigation_filename.c_str());
    std::remove(binningtable_filename.c_str());
    std::remove(config_filename.c_str());
}
