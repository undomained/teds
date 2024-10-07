# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# List of tests without the extension
set(tests unit)

# Create test targets
foreach(_test ${tests})
  add_executable(l1b_${_test} EXCLUDE_FROM_ALL ${_test}.cpp)
  target_link_libraries(l1b_${_test} PRIVATE tango_l1b Catch2::Catch2WithMain)
  target_compile_definitions(l1b_${_test} PRIVATE FIXTURE_DIR="${FIXTURE_DIR}")
  add_test(l1b_${_test} l1b_${_test})
  add_dependencies(test l1b_${_test})
endforeach()
