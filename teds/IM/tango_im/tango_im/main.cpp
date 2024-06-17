// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Read in settings, check for errors and inconsistencies in settings,
// and run the instrument model.
//
// When run without an argument print a verbose list of all parameters
// available to the user user. Otherwise run the instrument model with
// settings read from a YAML configuration file.

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
        tango::driver(settings, argc, argv);
    }
    return 0;
}
