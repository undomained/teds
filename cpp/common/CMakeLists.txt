# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Has to be shared in order to interface with Python
add_library(common SHARED
  # IO
  io.cpp
  time.cpp
  timer.cpp
  yaml.cpp
  settings.cpp
  # Math
  algorithm.cpp
  linalg.cpp
  linear_spline.cpp
  cubic_spline.cpp
  b_spline.cpp
  b_spline_2d.cpp
  fourier.cpp
  # Instrument
  ckd.cpp
  isrf.cpp
  binning_table.cpp
  # Geometry
  quaternion.cpp
  geometry.cpp
  dem.cpp
  solar_model.cpp)

# Process user defined include paths and libraries
target_include_directories(common PUBLIC ${INCLUDE_PATHS})
# Convert variables to list first
foreach(_var LIBRARIES LIBRARY_PATHS INCLUDE_PATHS)
  # If the variable was already given as list, do nothing
  list(LENGTH ${_var} _length)
  if (_length EQUAL 1)
    string(REPLACE " " ";" ${_var} ${${_var}})
  endif()
endforeach()
# Go through directories listed in LIBRARY_PATHS and turn entries of
# LIBRARIES into targets.
foreach(_lib ${LIBRARIES})
  # If target already exists skip it
  if (NOT TARGET ${_lib})
    find_library(LIB_FULLPATH ${_lib} HINTS ${LIBRARY_PATHS})
    if (LIB_FULLPATH)
      message(STATUS "Found ${LIB_FULLPATH}")
      add_library(${_lib} UNKNOWN IMPORTED)
      set_target_properties(${_lib} PROPERTIES
        IMPORTED_LOCATION ${LIB_FULLPATH})
      unset(LIB_FULLPATH CACHE)
    else()
      message(FATAL_ERROR "Could not find ${_lib}")
    endif()
  endif()
endforeach()
target_link_libraries(common PUBLIC ${LIBRARIES})

# C++20 features
target_compile_features(common PUBLIC cxx_std_20)

# Dependencies
include(FetchContent)
# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(common PUBLIC OpenMP::OpenMP_CXX)
# Python
find_package(Python3 COMPONENTS Interpreter Development NumPy)
if (TARGET Python3::NumPy)
  target_link_libraries(common PUBLIC Python3::Python Python3::NumPy)
endif()
# Logging
find_package(spdlog HINTS ${SPDLOG_PATH})
if (TARGET spdlog::spdlog)
  get_target_property(_loc spdlog::spdlog IMPORTED_LOCATION_NONE)
  message(STATUS "Found ${_loc}")
else()
  if (first_run)
    message("spdlog not found, downloading and building")
  endif()
  FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog
    GIT_TAG bae78f7b6c32fbfb97dd5751591734c12103fc1a)
  FetchContent_MakeAvailable(spdlog)
  install(DIRECTORY ${spdlog_SOURCE_DIR}/include/spdlog DESTINATION include)
  install(TARGETS spdlog EXPORT commonConfig)
endif()
target_link_libraries(common PUBLIC spdlog::spdlog)
# YAML
find_package(yaml-cpp 0.8 QUIET HINTS ${YAML_CPP_PATH})
if (TARGET yaml-cpp::yaml-cpp OR TARGET yaml-cpp)
  if (TARGET yaml-cpp::yaml-cpp)
    get_target_property(_loc yaml-cpp::yaml-cpp IMPORTED_LOCATION_RELEASE)
    target_link_libraries(common PUBLIC yaml-cpp::yaml-cpp)
  else()
    get_target_property(_loc yaml-cpp IMPORTED_LOCATION_RELEASE)
    target_link_libraries(common PUBLIC yaml-cpp)
  endif()
  message(STATUS "Found ${_loc}")
else()
  if (first_run)
    message("yaml-cpp not found, downloading and building")
  endif()
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
    GIT_TAG 28f93bdec6387d42332220afa9558060c8016795)
  FetchContent_MakeAvailable(yaml-cpp)
  target_link_libraries(common PUBLIC yaml-cpp)
  install(TARGETS yaml-cpp EXPORT commonConfig)
endif()
# Linear algebra library
find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(common PUBLIC Eigen3::Eigen)
# Fast Fourier transform library
find_package(pocketfft QUIET HINTS ${POCKETFFT_PATH})
if (TARGET pocketfft::pocketfft)
  get_target_property(_loc pocketfft::pocketfft INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found ${_loc}")
  target_link_libraries(common PUBLIC pocketfft::pocketfft)
else()
  if (first_run)
    message("pocketfft not found, downloading")
  endif()
  FetchContent_Declare(pocketfft
    GIT_REPOSITORY https://bitbucket.org/sron_earth/pocketfft
    GIT_TAG 35ef5eb45223b59c94146691a0ee50b8ec301909)
  FetchContent_MakeAvailable(pocketfft)
  get_target_property(_loc pocketfft INTERFACE_INCLUDE_DIRECTORIES)
  target_link_libraries(common PUBLIC pocketfft)
  install(TARGETS pocketfft EXPORT commonConfig)
endif()

# Retrieve paths of all linked libraries
get_target_property(_libraries common LINK_LIBRARIES)
set(library_paths "")
foreach(_lib ${_libraries})
  # Test if any of the following properties are defined. Those are all
  # candidates for the where the library is located.
  get_target_property(_loc ${_lib} IMPORTED_LOCATION)
  if (_loc MATCHES ".*-NOTFOUND")
    get_target_property(_loc ${_lib} IMPORTED_LOCATION_RELEASE)
  endif()
  if (_loc MATCHES ".*-NOTFOUND")
    get_target_property(_loc ${_lib} INTERFACE_INCLUDE_DIRECTORIES)
  endif()
  # If non of the above was found simply set the location to the
  # library name (basically remove the NOTFOUND part).
  if (_loc MATCHES ".*-NOTFOUND")
    set(_loc ${_lib})
  endif()
  # Clean the path if from a build instead of an install directory
  string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+).+" "\\1" _loc "${_loc}")
  list(APPEND library_paths ${_loc})
endforeach()
list(JOIN library_paths " " library_paths)
# Gather all compiler flags
set(all_flags ${CMAKE_CXX_FLAGS})
if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} buildtype)
  set(all_flags "${CMAKE_CXX_FLAGS_${buildtype}} ${all_flags}")
  string(STRIP ${all_flags} all_flags)
endif()
# Information about the environment shown during runtime
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GetGitRevisionDescription)
get_git_head_revision(_ GIT_COMMIT)
string(SUBSTRING ${GIT_COMMIT} 0 8 GIT_COMMIT_ABBREV)
target_compile_definitions(common PUBLIC
  TANGO_PROJECT_VERSION="${TANGO_VERSION}"
  TANGO_GIT_COMMIT_ABBREV="${GIT_COMMIT_ABBREV}"
  TANGO_CMAKE_HOST_SYSTEM="${CMAKE_HOST_SYSTEM}"
  TANGO_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
  TANGO_CXX_COMPILER_FLAGS="${all_flags}"
  TANGO_LIBRARIES="${library_paths}")

# Install
target_include_directories(common INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include>)
set(_headers ckd.h l1.h)
set_target_properties(common PROPERTIES PUBLIC_HEADER "${_headers}")
install(TARGETS common
  EXPORT commonConfig
  PUBLIC_HEADER DESTINATION include/common)
install(EXPORT commonConfig
  NAMESPACE common::
  DESTINATION lib/cmake/common)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/commonConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/commonConfigVersion.cmake
  DESTINATION lib/cmake/common)
export(EXPORT commonConfig NAMESPACE common::
  FILE ${PROJECT_BINARY_DIR}/commonConfig.cmake)
