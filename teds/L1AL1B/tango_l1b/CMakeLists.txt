# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# CMake
cmake_minimum_required(VERSION 3.25)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(Tools)

# Project
file(STRINGS ../../../project_version.txt TANGO_VERSION REGEX "^[0-9.]+")
project(tango_l1b VERSION ${TANGO_VERSION} LANGUAGES CXX)

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
option(BUILD_SHARED_LIBS
  "If enabled, builds a shared library, else static" OFF)

# Define main target
add_library(tango_l1b "")
add_executable(tango_l1b.x "")
target_link_libraries(tango_l1b.x PRIVATE tango_l1b)
include(GetGitRevisionDescription)
get_git_head_revision(_ GIT_COMMIT)
string(SUBSTRING ${GIT_COMMIT} 0 8 GIT_COMMIT_ABBREV)

# Process user defined include paths and libraries
target_include_directories(tango_l1b.x PUBLIC ${INCLUDE_PATHS})
# Convert variables to list first
foreach(_var LIBRARIES LIBRARY_PATHS INCLUDE_PATHS)
  convert_to_list(${_var})
endforeach()
generate_library_targets(LIBRARY_PATHS LIBRARIES)
target_link_libraries(tango_l1b.x PRIVATE ${LIBRARIES})

# Sources
add_subdirectory(tango_l1b)

# Retrieve locations from the list of supplied libraries and convert
# into a string.
get_target_property(all_libs tango_l1b LINK_LIBRARIES)
foreach(_lib ${LIBRARIES})
  get_target_property(_loc ${_lib} IMPORTED_LOCATION)
  # Clean the location if from a build instead of an install directory
  string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+).+" "\\1" _loc "${_loc}")
  list(APPEND all_libs ${_loc})
endforeach()
list(JOIN all_libs " " all_libs)
# Gather all compiler flags
get_all_flags(CXX all_flags)
target_compile_definitions(tango_l1b PUBLIC
  TANGO_PROJECT_VERSION="${TANGO_VERSION}"
  TANGO_GIT_COMMIT_ABBREV="${GIT_COMMIT_ABBREV}"
  TANGO_CMAKE_HOST_SYSTEM="${CMAKE_HOST_SYSTEM}"
  TANGO_EXECUTABLE="${CMAKE_CURRENT_BINARY_DIR}/tango_l1b"
  TANGO_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
  TANGO_CXX_COMPILER_FLAGS="${all_flags}"
  TANGO_LIBRARIES="${all_libs}")

# Print info
message(
  "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
get_all_flags(CXX all_flags)
message("Compiler flags: ${all_flags}")
message("Build configuration: ${CMAKE_BUILD_TYPE}")
message("Install location: ${CMAKE_INSTALL_PREFIX}")
