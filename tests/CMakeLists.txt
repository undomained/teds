# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

enable_testing()
check_include_coverage()

# Find and link against Catch2
set(CATCH2_PATH "" CACHE STRING "Path to Catch2 installation directory")
find_package(Catch2 3.3.2 QUIET HINTS ${CATCH2_PATH})
if (Catch2_FOUND)
  get_target_property(_loc Catch2::Catch2 INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found ${_loc}")
else()
  if (first_run)
    message("Catch2 >= 3.3.2 not found, downloading")
  endif()
  FetchContent_Declare(catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2 GIT_TAG v3.3.2)
  FetchContent_MakeAvailable(catch2)
endif()
add_custom_target(test COMMAND ctest --output-on-failure)

# Create target for code coverage if possible
if (INCLUDE_COVERAGE)
  find_program(_lcov lcov)
  if (_lcov)
    add_custom_target(coverage
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      COMMAND lcov -c -d . -o tracefile.info
      COMMAND
      perl -ne 'print if s|^SF:\(?!${PROJECT_SOURCE_DIR}\)|\\1|' tracefile.info
      | xargs lcov --ignore-errors unused,unused -r tracefile.info
      -o tracefile_filtered.info
      COMMAND genhtml tracefile_filtered.info -o out
      COMMAND cmake -E echo
      "Open ${PROJECT_BINARY_DIR}/out/index.html to see the detailed report"
      DEPENDS test)
  else()
    add_custom_target(coverage
      COMMAND cmake -E echo "Error: Cannot run coverage. Install lcov first.")
  endif()
else()
  add_custom_target(coverage
    COMMAND cmake -E echo
    "Error: Cannot run coverage. Enable INCLUDE_COVERAGE first.")
endif()

# Test data will be read from this directory by tests
set(FIXTURE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/common)

add_subdirectory(l1al1b/c++)
add_subdirectory(im/c++)
