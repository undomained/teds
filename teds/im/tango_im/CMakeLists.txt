# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

target_sources(tango_im PRIVATE
  driver.cpp
  driver_nitro.cpp
  settings_im.cpp
  uncalibration.cpp
)

file(GLOB_RECURSE SOURCES "*.cpp")
target_sources(tango_im.x PRIVATE main.cpp)


target_include_directories(tango_im INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
