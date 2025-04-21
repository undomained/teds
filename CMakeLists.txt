# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Top-level CMakeLists.txt for C++ components. L1B and L2 form the
# operational processor. The C++ targets can be built independently of
# the Python code, meaning this CMake project can be configured and
# built on its own.

# CMake
cmake_minimum_required(VERSION 3.22)

file(STRINGS pyproject.toml _version REGEX "^version *= *'[0-9.]+'$")
string(REGEX MATCH "[0-9.]+" TANGO_VERSION ${_version})
project(teds VERSION ${TANGO_VERSION} LANGUAGES CXX)

# The first_run variable is useful when attempting to detect the
# default value of a variable, but only during the first configuring
# and only if not already specified by the user (e.g., in the initial
# cache file).
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  set(first_run TRUE)
else()
  set(first_run FALSE)
endif()

# Options and variables
set(INCLUDE_PATHS "" CACHE STRING
  "Directories containing necessary header files")
set(LIBRARY_PATHS "" CACHE STRING
  "Directories containing libraries to be linked against")
set(LIBRARIES "" CACHE STRING "Libraries to be linked against")
set(SPDLOG_PATH "" CACHE STRING "Path to spdlog installation directory")
set(YAML_CPP_PATH "" CACHE STRING "Path to yaml-cpp installation directory")
set(POCKETFFT_PATH "" CACHE STRING "Path to pocketfft.h")
option(INCLUDE_COVERAGE "Define new target called 'coverage'" OFF)
set(COVERAGE_FLAGS "--coverage" CACHE STRING
  "Additional compiler flags used when running code coverage")
if (INCLUDE_COVERAGE)
  string(APPEND CMAKE_CXX_FLAGS " ${COVERAGE_FLAGS}")
endif()

# Libraries and executables to be built - one core library with common
# functionality, one library for each module (IM, L1B, L2), and one
# for tests.
add_subdirectory(cpp/common)
add_subdirectory(cpp/l1al1b)
add_subdirectory(cpp/im)
add_subdirectory(cpp/bindings)
add_subdirectory(tests)
