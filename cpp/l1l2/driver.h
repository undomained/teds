// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

namespace tango {

class SettingsL2;

// argc and argv are for generating the history attribute
auto driver(const SettingsL2& settings,
            const int argc = 0,
            const char* const argv[] = nullptr) -> void;

} // namespace tango
