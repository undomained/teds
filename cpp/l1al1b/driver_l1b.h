// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Entry point (driver function) for the Tango L1A-L1B calibration
// software. It initializes the CKD and runs the processing steps.

#pragma once

namespace tango {

class SettingsL1B;

// argc and argv are for generating the NetCDF history attribute
auto driverL1B(const SettingsL1B& settings,
               const int argc = 0,
               const char* const argv[] = nullptr) -> void;

} // namespace tango
