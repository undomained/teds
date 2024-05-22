# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Add source files and definitions (e.g. project version) to the
# tango_l1b target.

target_sources(tango_l1b PRIVATE
  driver.cpp
  processor.cpp
  array.cpp
  bspline.cpp
  fourier.cpp
  functions.cpp
  inversion.cpp
  lininv.cpp
  matrix.cpp
  netcdf_object.cpp
  planet.cpp
  settings.cpp
  vector.cpp
  logger.cpp
  batch.cpp
  binningtable.cpp
  ckd.cpp
  l1a.cpp
  l1a_file_metadata.cpp
  l1a_flight.cpp
  l1a_manager.cpp
  l1b.cpp
  l1x.cpp
  settings_main.cpp
  settings_proc.cpp
  cubic_spline.cpp
  parallel.cpp
)
target_sources(tango_l1b.x PRIVATE main.cpp)

# C++20 features
target_compile_features(tango_l1b PUBLIC cxx_std_20)

# Exporting
target_include_directories(tango_l1b INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include>)
set(_headers
  ckh.h)
set_target_properties(tango_l1b PROPERTIES PUBLIC_HEADER "${_headers}")
install(TARGETS tango_l1b
  EXPORT tango_l1bConfig
  PUBLIC_HEADER DESTINATION include/tango_l1b)
install(EXPORT tango_l1bConfig
  NAMESPACE tango_l1b::
  DESTINATION lib/cmake/tango_l1b)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/tango_l1bConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/tango_l1bConfigVersion.cmake
  DESTINATION lib/cmake/tango_l1b)
export(EXPORT tango_l1bConfig NAMESPACE tango_l1b::
  FILE ${PROJECT_BINARY_DIR}/tango_l1bConfig.cmake)
