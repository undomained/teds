// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Timer class for this project. Only works on the master OpenMP
// thread.

#pragma once

#include <chrono>

namespace tango {

class Timer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> wall_timestamp;
    double total_wall_time {};

public:
    Timer() = default;
    auto start() -> void;
    auto stop() -> void;
    // Return total wall time
    [[nodiscard]] auto time() const -> double;
};

} // namespace tango
