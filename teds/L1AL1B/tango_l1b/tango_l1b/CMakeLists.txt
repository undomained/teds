# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

# Add source files and definitions (e.g. project version) to the
# tango_l1b.x target.

target_sources(tango_l1b.x PRIVATE
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
  main.cpp
  settings_main.cpp
  settings_proc.cpp
  cubic_spline.cpp
  parallel.cpp
)

target_include_directories(tango_l1b.x PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR})

# Link against a linear algebra library
if (LINALG_FLAVOR MATCHES "mkl")
  # If LINALG_FLAVOR has been set then use that library
  find_package(MKL REQUIRED)
  target_include_directories(tango_l1b.x PUBLIC ${MKL_INCLUDE_DIR})
  target_link_libraries(tango_l1b.x PUBLIC ${MKL_LIBRARIES})
  target_compile_definitions(tango_l1b.x PUBLIC USE_MKL)
  set(linalg_library_name "MKL")
else()
  # Otherwise attempt to detect if MKL is used
  try_compile(LINALG_TEST ${CMAKE_CURRENT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/tests/MKL.cpp
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${INCLUDE_PATHS}"
    "-DLINK_DIRECTORIES=${LIBRARY_PATHS}"
    "-DLINK_LIBRARIES=${LIBRARIES}")
  if (LINALG_TEST)
    target_compile_definitions(tango_l1b.x PUBLIC USE_MKL)
    set(linalg_library_name "MKL")
  else()
    set(linalg_library_name "system default")
  endif()
endif()
message("Using linear algebra library: ${linalg_library_name}")
