// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Time and date related functions

#pragma once

#include <string>

namespace tango {

auto getDate() -> std::string;

// Format YYYY-mm-ddTHH:MM:SSZ
auto getDateAndTime() -> std::string;

} // namespace tango
