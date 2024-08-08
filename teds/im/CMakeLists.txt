# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# CMake
cmake_minimum_required(VERSION 3.22)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/../l1al1b/cmake)
include(Tools)

# Project
file(STRINGS ../../project_version.txt TANGO_VERSION REGEX "^[0-9.]+")
project(tango_im VERSION ${TANGO_VERSION} LANGUAGES CXX)

# Options and variables
set(INCLUDE_PATHS "" CACHE STRING
  "Directories containing necessary header files")
set(LIBRARY_PATHS "" CACHE STRING
  "Directories containing libraries to be linked against")
set(LIBRARIES "" CACHE STRING "Libraries to be linked against")
set(TANGO_L1B_PATH "" CACHE STRING "Path to the Tango L1B library")
option(BUILD_SHARED_LIBS
  "If enabled, dependencies are built as shared libraries by default" OFF)

# Define main target
add_executable(tango_im.x "")
target_compile_features(tango_im.x PUBLIC cxx_std_20)
include(GetGitRevisionDescription)
get_git_head_revision(_ GIT_COMMIT)
string(SUBSTRING ${GIT_COMMIT} 0 8 GIT_COMMIT_ABBREV)

# Sources
add_subdirectory(tango_im)

# Dependencies
find_package(OpenMP REQUIRED)
find_package(Threads)
find_package(spdlog QUIET)
find_package(yaml-cpp 0.7 QUIET)
find_package(tango_l1b REQUIRED HINTS ${TANGO_L1B_PATH})
target_link_libraries(tango_im.x PRIVATE tango_l1b::tango_l1b)

# Process user defined include paths and libraries
target_include_directories(tango_im.x PUBLIC ${INCLUDE_PATHS})
# Convert variables to list first
foreach(_var LIBRARIES LIBRARY_PATHS INCLUDE_PATHS)
  convert_to_list(${_var})
endforeach()
generate_library_targets(LIBRARY_PATHS LIBRARIES)
target_link_libraries(tango_im.x PRIVATE ${LIBRARIES})

# Retrieve locations from the list of supplied libraries and convert
# into a string.
get_target_property(all_libs tango_im.x LINK_LIBRARIES)
foreach(_lib ${LIBRARIES})
  get_target_property(_loc ${_lib} IMPORTED_LOCATION)
  # Clean the location if from a build instead of an install directory
  string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+).+" "\\1" _loc "${_loc}")
  list(APPEND all_libs ${_loc})
endforeach()
list(JOIN all_libs " " all_libs)
# Gather all compiler flags
get_all_flags(CXX all_flags)
target_compile_definitions(tango_im.x PRIVATE
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
