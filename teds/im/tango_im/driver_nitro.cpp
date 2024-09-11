// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver_nitro.h"

#include "settings_im.h"
#include "uncalibration.h"

#include <netcdf>
#include <spdlog/spdlog.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/io.h>
#include <tango_l1b/l1.h>
#include <tango_l1b/l1_measurement.h>
#include <tango_l1b/timer.h>
#include <tango_l1b/algorithms/build_algo.h>
#include <tango_l1b/algorithms/base_algo.h>

namespace tango {


// Meta data such as the exposure time are read from the L1A input
// file but ignored by the instrument model. Instead, they are set by
// user settings. The L1A-L1B processor, however, only reads from the
// input file because after the IM run these parameter should be
// fixed.
static auto setL1Meta(const SettingsIM& settings,
                      std::vector<L1>& l1_products) -> void
{
    for (auto& l1 : l1_products) {
        l1.binning_table_id =
          static_cast<uint8_t>(settings.detector.binning_table_id);
        l1.nr_coadditions =
          static_cast<uint16_t>(settings.detector.nr_coadditions);
        l1.exposure_time = settings.detector.exposure_time;
    }
}

auto driver_nitro(const SettingsIM& settings,
            const int argc,
            const char* const argv[]) -> void
{
    // Set up loggers and print general information
    initLogging(false);
    printHeading("Tango instrument model", false);
    printSystemInfo(TANGO_PROJECT_VERSION,
                    TANGO_GIT_COMMIT_ABBREV,
                    TANGO_CMAKE_HOST_SYSTEM,
                    TANGO_EXECUTABLE,
                    TANGO_CXX_COMPILER,
                    TANGO_CXX_COMPILER_FLAGS,
                    TANGO_LIBRARIES,
                    settings.io.binning_table);

    // Read in the CKD
    printHeading("Reading CKD and input data");
    const CKD ckd { settings.io.ckd };

    printHeading("Initialize the binning table");
    // Initialize the binning table
    const BinningTable binning_table { ckd.n_detector_rows,
                                       ckd.n_detector_cols,
                                       settings.io.binning_table,
                                       settings.detector.binning_table_id };

    const std::string& config = settings.getConfig();

    L1Measurement l1_measurement(settings.io.sgm, settings.image_start, settings.image_end.value_or(fill::i), config);

    const std::string& proctable_file = settings.proctable.file;
    spdlog::info("proctable_file: {} ", proctable_file);

    const std::string& algo_list_name = settings.proctable.algo_list;
    spdlog::info("algo_list_name: {} ", algo_list_name);

    YAML::Node proctable = YAML::LoadFile(proctable_file);
    YAML::Node algo_list = proctable[algo_list_name];

    // Run the forward model (main loop)
    printHeading("Instrument Model");
    std::array<Timer, static_cast<int>(ProcLevel::n_levels)> timers {};
    Timer timer_total {};
    timer_total.start();

//    #pragma omp parallel for schedule(dynamic)

    for (int i_alt = 0; i_alt < static_cast<int>(l1_measurement.size()); ++i_alt) {
        printPercentage(i_alt, l1_measurement.size(), "Processing scenes");

        auto& l1 { l1_measurement[i_alt] };
        
        // Initialize pixel mask
        l1.pixel_mask = ckd.pixel_mask;
        
        std::string i_alt_msg = " Processing image [" + std::to_string(i_alt) +  "] ";
        spdlog::info("{:=^30}", i_alt_msg);
        int i_algo = 0;
        for (YAML::const_iterator it=algo_list.begin(); it!=algo_list.end();it++){
            BuildAlgo algo_builder;
            BaseAlgo* algo = algo_builder.CreateAlgo(it->as<std::string>());
            
            // Inverse operators
            algo->setModelType("IM");

            if (algo) {
                timers[static_cast<int>(i_algo)].start();
                if (algo->algoCheckInput(ckd, l1)) {
                    spdlog::info("{: ^30}", algo->getName()); // Remove this later
                    algo->algoExecute(ckd, l1);
                }
                timers[static_cast<int>(i_algo)].stop();
            }
            i_algo+=1;    
        }
    }
    timer_total.stop();
    spdlog::info("Processing images 100.0%");

    // For writing to output switch out the LBL wavelength grid
    //*l1_products.front().wavelength = ckd.wave.wavelength;

    // Overview of timings
    spdlog::info("");
    spdlog::info("Timings:");
    int i_algo = 0;
    for (YAML::const_iterator it=algo_list.begin(); it!=algo_list.end();it++){
        std::string my_name = it->as<std::string>();
        spdlog::info("   {:20}: {:8.3f} s",
                    my_name, timers[static_cast<int>(i_algo)].time());
        i_algo += 1;
    }
    spdlog::info("   {:20}: {:8.3f} s", "Total",timer_total.time());

    printHeading("Writing file");

    std::string level = "L1A";
    if (algo_list.size() <= 1){
        level = "SGM";
    }
    l1_measurement.write(settings.io.l1a, level, settings.getConfig(), argc, argv);

    printHeading("Success");
}

} // namespace tango
