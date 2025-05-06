# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

add_library(tango_im SHARED
  driver_im.cpp
  settings_im.cpp
  forward_models.cpp)
target_link_libraries(tango_im PUBLIC common)

add_executable(tango_im.x main.cpp)
target_link_libraries(tango_im.x PRIVATE tango_im)

target_compile_definitions(tango_im PRIVATE
  TANGO_EXECUTABLE="${CMAKE_CURRENT_BINARY_DIR}/tango_im.x")

# Install
target_include_directories(tango_im INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)
