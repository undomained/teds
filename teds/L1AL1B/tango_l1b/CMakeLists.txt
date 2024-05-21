# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# CMake
cmake_minimum_required(VERSION 3.25)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include(Tools)

# When using newer version of CMake some commands behave
# differently. This conditional helps to guard against that.
if (CMAKE_VERSION VERSION_LESS 3.13.5)
  include(redefine_target_sources)
else()
  cmake_policy(VERSION 3.13.5)
endif()

# Project
project(tango_l1b VERSION 1.0 LANGUAGES CXX)

# Options and variables
set(INCLUDE_PATHS "" CACHE STRING
  "Directories containing necessary header files")
set(LIBRARY_PATHS "" CACHE STRING
  "Directories containing libraries to be linked against")
set(LIBRARIES "" CACHE STRING "Libraries to be linked against")
set(LINALG_FLAVOR "" CACHE STRING "Choose the linear algebra library")
set(SPDLOG_PATH "" CACHE STRING "Path to spdlog installation directory")
set(GADFIT_PATH "" CACHE STRING "Path to GADfit installation directory")

# Attempt to detect the default values of some variables, but only
# during the first configuring and only if not already specified by
# the user (e.g., in the initial cache file).
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  # Detect USE_MPI
  try_compile(MPI_OK ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/tests/MPI.cpp)
  option(USE_MPI "Whether to compile with MPI support" ${MPI_OK})
endif()

# Define main target
add_executable(tango_l1b.x "")
target_compile_features(tango_l1b.x PUBLIC cxx_std_20)

# MPI (Test program does not need MPI, maybe later)
if (USE_MPI)
  target_compile_definitions(tango_l1b.x PUBLIC USE_MPI)
  message("Using MPI")
else()
  message("Not using MPI")
endif()

# Process user defined include paths and libraries
target_include_directories(tango_l1b.x PUBLIC ${INCLUDE_PATHS})
# Convert variables to list first
foreach(_var LIBRARIES LIBRARY_PATHS INCLUDE_PATHS)
  convert_to_list(${_var})
endforeach()
generate_library_targets(LIBRARY_PATHS LIBRARIES)
target_link_libraries(tango_l1b.x PRIVATE ${LIBRARIES})

# Print info
message(
  "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
get_all_flags(CXX all_flags)
message("Compiler flags: ${all_flags}")
message("Build configuration: ${CMAKE_BUILD_TYPE}")
message("Install location: ${CMAKE_INSTALL_PREFIX}")

# Sources
add_subdirectory(tango_l1b)
add_subdirectory(external/gadfit)
target_link_libraries(tango_l1b.x PUBLIC gadfit)
