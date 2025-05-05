# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

add_library(tango_l2 SHARED
  driver.cpp
  settings_l2.cpp
  atmosphere.cpp
  optic_abs_prop.cpp
  read_sun_spectrum.cpp
  gauss_newton.cpp
  surface.cpp
  nonscat_fwd_model.cpp
)
target_link_libraries(tango_l2 PUBLIC common)

add_executable(tango_l2.x main.cpp)
target_link_libraries(tango_l2.x PRIVATE tango_l2)

target_compile_definitions(tango_l2 PRIVATE
  TANGO_EXECUTABLE="${CMAKE_CURRENT_BINARY_DIR}/tango_l2.x")

# Install
target_include_directories(tango_l2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
