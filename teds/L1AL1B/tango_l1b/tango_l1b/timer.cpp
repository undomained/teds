// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "timer.h"

#include <omp.h>

namespace tango {

auto Timer::start() -> void
{
    if (omp_get_thread_num() == 0) {
        wall_timestamp = std::chrono::high_resolution_clock::now();
    }
}

auto Timer::stop() -> void
{
    if (omp_get_thread_num() == 0) {
        total_wall_time +=
          std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::high_resolution_clock::now() - wall_timestamp)
            .count();
    }
}

[[nodiscard]] auto Timer::time() const -> double
{
    return total_wall_time;
}

} // namespace tango
