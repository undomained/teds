// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "driver.h"
#include "settings_l2.h"

#include <iostream>

auto main(int argc, char* argv[]) -> int
{
    std::cout.precision(16);
    if (argc == 1) {
        // When called without an argument print the default
        // configuration.
        std::cout << "%YAML 1.2\n---\n" << tango::SettingsL2 {}.c_str() << '\n';
    } else {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        tango::SettingsL2 settings { argv[1] };
        settings.init();
        tango::driver(settings, argc, argv);
    }
    return 0;
}
