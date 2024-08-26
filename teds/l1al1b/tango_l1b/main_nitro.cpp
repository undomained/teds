// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Read in settings, check for errors and inconsistencies in settings,
// and run the L1A-L1B processor.
//
// When run without an argument print a verbose list of all parameters
// available to the user user. Otherwise run the processor with
// settings read from a YAML configuration file.

#include "driver_nitro.h"
#include "settings_l1b.h"

#include <iostream>

auto main(int argc, char* argv[]) -> int
{
    if (argc == 1) {
        // When called without an argument print the default
        // configuration.
        std::cout << "%YAML 1.2\n---\n"
                  << tango::SettingsL1B {}.c_str() << '\n';
    } else {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        tango::SettingsL1B settings { argv[1] };
        settings.init();
        driver(settings, argc, argv);
    }
    return 0;
}
