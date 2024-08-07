// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "time.h"

#include <ctime>

namespace tango {

constexpr int date_size { 100 };

auto getDate() -> std::string
{
    std::time_t t { std::time(nullptr) };
    char date[date_size] {};
    std::strftime(
      date, date_size * sizeof(char), "%Y %B %d %a UTC%z", std::localtime(&t));
    return date;
}

auto getDateAndTime() -> std::string
{
    std::time_t t { std::time(nullptr) };
    char date_and_time[date_size] {};
    std::strftime(
      date_and_time, date_size * sizeof(char), "%FT%TZ", std::gmtime(&t));
    return date_and_time;
}

} // namespace tango
