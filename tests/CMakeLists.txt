# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

enable_testing()

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
      COMMAND lcov --directory . --capture --output-file coverage.info
      COMMAND genhtml --ignore-errors inconsistent --demangle-cpp -o coverage coverage.info
      COMMAND cmake -E echo
      "Open ${PROJECT_BINARY_DIR}/coverage/index.html and click on tango_im and tango_l1b entries to see the detailed reports"
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

add_subdirectory(cpp/im)
add_subdirectory(cpp/l1al1b)
