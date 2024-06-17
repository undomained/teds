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
set(SPDLOG_PATH "" CACHE STRING "Path to spdlog installation directory")
set(YAML_CPP_PATH "" CACHE STRING "Path to yaml-cpp installation directory")
set(POCKETFFT_PATH "" CACHE STRING "Path to pocketfft.h")
option(BUILD_SHARED_LIBS
  "If enabled, builds a shared library, else static" OFF)
option(INCLUDE_COVERAGE "Defines new target called 'coverage'" OFF)
set(COVERAGE_FLAGS "--coverage" CACHE STRING
  "Additional compiler flags used when running code coverage")
if (INCLUDE_COVERAGE)
  string(APPEND CMAKE_CXX_FLAGS " ${COVERAGE_FLAGS}")
endif()

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

# Dependencies
include(FetchContent)
find_package(OpenMP REQUIRED)
target_link_libraries(tango_l1b PUBLIC OpenMP::OpenMP_CXX)
find_package(spdlog HINTS ${SPDLOG_PATH})
if (TARGET spdlog::spdlog)
  message("target found")
else()
  if (first_run)
    message("spdlog not found, downloading")
  endif()
  FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog
    GIT_TAG bae78f7b6c32fbfb97dd5751591734c12103fc1a)
  FetchContent_MakeAvailable(spdlog)
  install(DIRECTORY ${spdlog_SOURCE_DIR}/include/spdlog DESTINATION include)
  install(TARGETS spdlog EXPORT tango_l1bConfig)
endif()
target_link_libraries(tango_l1b PUBLIC spdlog::spdlog)
find_package(yaml-cpp 0.7 QUIET HINTS ${YAML_CPP_PATH})
if (TARGET yaml-cpp::yaml-cpp OR TARGET yaml-cpp)
  if (TARGET yaml-cpp::yaml-cpp)
    get_target_property(_loc yaml-cpp::yaml-cpp IMPORTED_LOCATION_RELEASE)
  else()
    get_target_property(_loc yaml-cpp IMPORTED_LOCATION_RELEASE)
  endif()
  message(STATUS "Found ${_loc}")
  target_link_libraries(tango_l1b PUBLIC yaml-cpp)
else()
  if (first_run)
    message("yaml-cpp >= 0.7 not found, downloading")
  endif()
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
    GIT_TAG 5057a5329332993fa5618f4c7f1aff4c31c21413)
  FetchContent_MakeAvailable(yaml-cpp)
  target_link_libraries(tango_l1b PUBLIC yaml-cpp)
  install(TARGETS yaml-cpp EXPORT tango_l1bConfig)
endif()
find_package(pocketfft QUIET HINTS ${POCKETFFT_PATH})
if (TARGET pocketfft::pocketfft)
  get_target_property(_loc pocketfft::pocketfft INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found ${_loc}")
  target_link_libraries(tango_l1b PUBLIC pocketfft::pocketfft)
else()
  if (first_run)
    message("pocketfft not found, downloading")
  endif()
  FetchContent_Declare(pocketfft
    GIT_REPOSITORY https://bitbucket.org/sron_earth/pocketfft
    GIT_TAG 35ef5eb45223b59c94146691a0ee50b8ec301909)
  FetchContent_MakeAvailable(pocketfft)
  get_target_property(_loc pocketfft INTERFACE_INCLUDE_DIRECTORIES)
  target_link_libraries(tango_l1b PUBLIC pocketfft)
  install(TARGETS pocketfft EXPORT tango_l1bConfig)
endif()

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
target_compile_definitions(tango_l1b PRIVATE
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
