// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Entry point (driver function) for the Tango instrument model. It
// initializes the CKD and runs the steps to "uncalibrated" spectra or
// detector images.

#pragma once

namespace tango {

class SettingsIM;

// argc and argv are for generating the history attribute
auto driverIM(const SettingsIM& settings,
              const int argc = 0,
              const char* const argv[] = nullptr) -> void;

} // namespace tango
