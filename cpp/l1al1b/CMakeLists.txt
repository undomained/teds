# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

add_library(tango_l1b SHARED
  driver.cpp
  settings_l1b.cpp
  calibration.cpp)
target_link_libraries(tango_l1b PUBLIC common)

add_executable(tango_l1b.x main.cpp)
target_link_libraries(tango_l1b.x PRIVATE tango_l1b)

target_compile_definitions(tango_l1b PRIVATE
  TANGO_EXECUTABLE="${CMAKE_CURRENT_BINARY_DIR}/tango_l1b.x")

# Install
target_include_directories(tango_l1b INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
