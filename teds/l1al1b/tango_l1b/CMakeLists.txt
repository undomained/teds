# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

file(GLOB_RECURSE SOURCES "*.cpp")

target_sources(tango_l1b PRIVATE ${SOURCES})
target_sources(tango_l1b.x PRIVATE main.cpp)
#target_sources(tango_l1b.x PRIVATE main_nitro.cpp)

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
