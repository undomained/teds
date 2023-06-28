# Licensed under the Apache License, Version 2.0 (the "License"); you
# may not use this file except in compliance with the License.  You
# may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.  See the License for the specific language governing
# permissions and limitations under the License.

if (CMAKE_VERSION VERSION_LESS 3.13.5)
  include(redefine_target_sources)
else()
  cmake_policy(VERSION 3.13.5)
endif()

target_sources(gadfit PRIVATE
  automatic_differentiation.cpp
  fit_function.cpp
  lm_solver.cpp
  exceptions.cpp
  parallel.cpp
  numerical_integration.cpp
  timer.cpp
  )

# First run
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  # Automatically detect USE_MPI
  try_compile(MPI_OK ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/tests/MPI.cpp)
  option(USE_MPI "Whether to compile with MPI support" ${MPI_OK})
endif()

# Linear algebra library
if (LINALG_FLAVOR MATCHES "mkl")
  find_package(MKL REQUIRED)
  target_include_directories(gadfit PRIVATE ${MKL_INCLUDE_DIR})
  target_link_libraries(gadfit PRIVATE ${MKL_LIBRARIES})
  target_sources(gadfit PRIVATE cblas_wrapper_mkl.cpp)
  set(linalg_library_name "MKL")
elseif (LINALG_FLAVOR MATCHES "none")
  target_sources(gadfit PRIVATE cblas_wrapper_fallback.cpp)
  set(linalg_library_name "fallback")
else()
  try_compile(LINALG_TEST ${CMAKE_CURRENT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/tests/MKL.cpp
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${INCLUDE_PATHS}"
    "-DLINK_DIRECTORIES=${LIBRARY_PATHS}"
    "-DLINK_LIBRARIES=${LIBRARIES}")
  if (LINALG_TEST)
    target_sources(gadfit PRIVATE cblas_wrapper_mkl.cpp)
    set(linalg_library_name "MKL")
  else()
    target_sources(gadfit PRIVATE cblas_wrapper_fallback.cpp)
    set(linalg_library_name "fallback")
  endif()
endif()
message("Using linear algebra library: ${linalg_library_name}")

# MPI
if (USE_MPI)
  target_compile_definitions(gadfit PUBLIC USE_MPI)
  message("Using MPI")
else()
  message("Not using MPI")
endif()

# C++ features
if (NOT CMAKE_VERSION VERSION_LESS 3.8.2)
  target_compile_features(gadfit PUBLIC cxx_std_17)
endif()

# Exporting
target_include_directories(gadfit INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)
set(_headers
  lm_solver.h automatic_differentiation.h exceptions.h fit_function.h
  parallel.h numerical_integration.h timer.h)
set_target_properties(gadfit PROPERTIES PUBLIC_HEADER "${_headers}")
if (NOT CMAKE_VERSION VERSION_LESS 3.14.7)
  install(TARGETS gadfit
    EXPORT gadfitConfig
    PUBLIC_HEADER DESTINATION include/gadfit)
  install(EXPORT gadfitConfig NAMESPACE gadfit:: DESTINATION lib/cmake/gadfit)
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/gadfitConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
  install(FILES ${PROJECT_BINARY_DIR}/gadfitConfigVersion.cmake
    DESTINATION lib/cmake/gadfit)
  export(EXPORT gadfitConfig NAMESPACE gadfit::
    FILE ${PROJECT_BINARY_DIR}/gadfitConfig.cmake)
endif()
