# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# CMake
cmake_minimum_required(VERSION 3.25)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../L1AL1B/tango_l1b/cmake)
include(Tools)

# Project
file(STRINGS ../../../project_version.txt TANGO_VERSION REGEX "^[0-9.]+")
project(tango_im VERSION ${TANGO_VERSION} LANGUAGES CXX)

# Options and variables
set(INCLUDE_PATHS "" CACHE STRING
  "Directories containing necessary header files")
set(LIBRARY_PATHS "" CACHE STRING
  "Directories containing libraries to be linked against")
set(LIBRARIES "" CACHE STRING "Libraries to be linked against")
option(BUILD_SHARED_LIBS
  "If enabled, dependencies are built as shared libraries by default" OFF)

# Define main target
add_executable(tango_im.x "")
target_compile_features(tango_im.x PUBLIC cxx_std_20)

# Process user defined include paths and libraries
target_include_directories(tango_im.x PUBLIC ${INCLUDE_PATHS})
# Convert variables to list first
foreach(_var LIBRARIES LIBRARY_PATHS INCLUDE_PATHS)
  convert_to_list(${_var})
endforeach()
generate_library_targets(LIBRARY_PATHS LIBRARIES)
target_link_libraries(tango_im.x PRIVATE ${LIBRARIES})

# Sources
add_subdirectory(tango_im)

# Print info
message(
  "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
get_all_flags(CXX all_flags)
message("Compiler flags: ${all_flags}")
message("Build configuration: ${CMAKE_BUILD_TYPE}")
message("Install location: ${CMAKE_INSTALL_PREFIX}")
