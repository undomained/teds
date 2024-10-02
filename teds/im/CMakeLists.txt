# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Most of the build environment is defined in the top-level
# CMakeLists.txt of the L1A-L1B processor. The instrument model uses
# L1A-L1B as a dependency and only defines properties (flags,
# libraries, ..) that are not required for L1A-L1B.

# Define main target
add_library(tango_im "")
add_executable(tango_im.x "")
target_link_libraries(tango_im.x PRIVATE tango_im)

# Sources
add_subdirectory(tango_im)

# Dependencies
target_link_libraries(tango_im PUBLIC tango_l1b)

# Information about the environment that shown when running the
# instrument model
get_target_property(all_flags tango_l1b TEDS_ALL_FLAGS)
get_target_property(all_libs tango_l1b TEDS_ALL_LIBS)
get_target_property(commit_hash tango_l1b TEDS_COMMIT_HASH)
target_compile_definitions(tango_im PRIVATE
  TANGO_PROJECT_VERSION="${TANGO_VERSION}"
  TANGO_GIT_COMMIT_ABBREV="${commit_hash}"
  TANGO_CMAKE_HOST_SYSTEM="${CMAKE_HOST_SYSTEM}"
  TANGO_EXECUTABLE="${CMAKE_CURRENT_BINARY_DIR}/tango_l1b"
  TANGO_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
  TANGO_CXX_COMPILER_FLAGS="${all_flags}"
  TANGO_LIBRARIES="${all_libs}")
