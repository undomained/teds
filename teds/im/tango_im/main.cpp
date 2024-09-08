// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Read in settings, check for errors and inconsistencies in settings,
// and run the instrument model.
//
// When run without an argument print a verbose list of all parameters
// available to the user user. Otherwise run the instrument model with
// settings read from a YAML configuration file.

#include "driver_nitro.h"
#include "driver.h"
#include "settings_im.h"

#include <iostream>

auto main(int argc, char* argv[]) -> int
{
    if (argc == 1) {
        // When called without an argument print the default
        // configuration.
        std::cout << "%YAML 1.2\n---\n" << tango::SettingsIM {}.c_str() << '\n';
    } else {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        tango::SettingsIM settings { argv[1] };
        settings.init();
        const std::string instrument = settings.instrument;
        const std::string nitro = "nitro";
        if (instrument.compare(nitro) == 0){
            std::cout << "NITRO DRIVER" << std::endl;
            driver_nitro(settings, argc, argv);
        } else {
            std::cout << "CARBO DRIVER" << std::endl;
            driver(settings, argc, argv);
        }
    }
    return 0;
}
