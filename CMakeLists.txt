# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Top-level CMakeLists.txt for both the L1A-L1B processor, the
# instrument model, and tests. The two modules and the tests are in
# different directories but depend on each other. This wrapper brings
# them into a single CMake project which simplifies the build
# process.
#
# The L1A-L1B module is currently the only one that forms part of the
# operational processor. Its source directory contains full CMake
# support and it can be built independently as a CMake project. The
# other modules depend on the L1A-L1B processor and cannot be built
# independently.

# CMake
cmake_minimum_required(VERSION 3.22)

project(tango VERSION 0.1 LANGUAGES CXX)

add_subdirectory(teds/l1al1b)
add_subdirectory(teds/im)
add_subdirectory(tests)
