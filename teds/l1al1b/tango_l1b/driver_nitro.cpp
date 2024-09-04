// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver_nitro.h"

#include "algorithms/build_algo.h"
#include "algorithms/base_algo.h"
#include "binning_table.h"
#include "calibration.h"
#include "ckd.h"
#include "io.h"
#include "l1.h"
#include "settings_l1b.h"
#include "timer.h"
#include <yaml-cpp/yaml.h>

#include <spdlog/spdlog.h>

namespace tango {

auto driver_nitro(const SettingsL1B& settings,
            const int argc,
            const char* const argv[]) -> void
{
    // Set up loggers and print general information
    initLogging(false);
    printHeading("Tango L1B processor", false);
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
    CKD ckd(settings.io.ckd, settings.swath.spectrum_width);

    // Initialize L1 products by reading all L1A data (everything is
    // stored in memory).
    std::vector<L1> l1_products {};
    readL1(settings.io.l1a,
           settings.image_start,
           settings.image_end.value_or(fill::i),
           l1_products);

    // Initialize the binning table and bin the CKD
    const BinningTable binning_table {
        ckd.n_detector_rows,
        ckd.n_detector_cols,
        settings.io.binning_table,
        static_cast<int>(l1_products.front().binning_table_id)
    };
    if (settings.unbinning == Unbin::none) {
        binning_table.bin(ckd.pixel_mask);
        binning_table.bin(ckd.dark.offset);
        binning_table.bin(ckd.dark.current);
        binning_table.bin(ckd.noise.g);
        binning_table.bin(ckd.noise.n2);
        binning_table.bin(ckd.prnu.prnu);
        binning_table.binPixelIndices(ckd.swath.pix_indices);
    }

    // Run retrieval
    printHeading("Retrieval");
    Timer timer_total {};
    timer_total.start();

    const std::string& config = settings.getConfig();

    const std::string& proctable_file = settings.proctable.file;
    spdlog::info("proctable_file: {} ", proctable_file);

    const std::string& algo_list_name = settings.proctable.algo_list;
    spdlog::info("algo_list_name: {} ", algo_list_name);


    YAML::Node proctable = YAML::LoadFile(proctable_file);
    YAML::Node algo_list = proctable[algo_list_name];

    const int algo_list_length = algo_list.size();
// For some reason array is not happy with algo_list_length
//    std::array<Timer, algo_list_length> timers {};
    std::array<Timer, 20> timers {};

    #pragma omp parallel for schedule(dynamic)
    for (int i_alt = 0; i_alt < static_cast<int>(l1_products.size()); ++i_alt) {
        printPercentage(i_alt, l1_products.size(), "Processing images");
        auto& l1 { l1_products[i_alt] };
        
        // Initialize pixel mask
        l1.pixel_mask = ckd.pixel_mask;
        
        std::string i_alt_msg = " Processing image [" + std::to_string(i_alt) +  "] ";
        spdlog::info("{:=^30}", i_alt_msg);
        int i_algo = 0;
        for (YAML::const_iterator it=algo_list.begin(); it!=algo_list.end();it++){
            spdlog::info("{: ^30}", it->as<std::string>()); // Remove this later
            BuildAlgo algo_builder;
            BaseAlgo* algo = algo_builder.CreateAlgo(it->as<std::string>());
            if (algo) {
                timers[static_cast<int>(i_algo)].start();
                if (algo->algoCheckInput(ckd, l1)) {
                    algo->algoExecute(ckd, l1);
                }
                timers[static_cast<int>(i_algo)].stop();
            }
            i_algo+=1;    
        }
    }
    timer_total.stop();
    spdlog::info("Processing images 100.0%");

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

    writeL1product(settings.io.l1b, "L1B", settings.getConfig(), l1_products, argc, argv);

    printHeading("Success");
}

} // namespace tango
