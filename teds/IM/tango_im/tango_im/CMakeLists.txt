# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

target_sources(tango_im.x PRIVATE
  array.cpp
  binningtable.cpp
  bspline.cpp
  ckd.cpp
  fourier.cpp
  frame.cpp
  functions.cpp
  l1x.cpp
  l1x_inputfile.cpp
  logger.cpp
  main.cpp
  netcdf_object.cpp
  planet.cpp
  program.cpp
  settings.cpp
  settings_geo.cpp
  settings_isrf.cpp
  settings_l1b.cpp
  settings_main.cpp
  settings_noise.cpp
  utc.cpp
  vector.cpp
  cubic_spline.cpp
  linear_spline.cpp
)

get_target_property(_sources tango_im.x SOURCES)

target_include_directories(tango_im.x PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR})

# Link against a linear algebra library
if (LINALG_FLAVOR MATCHES "mkl")
  # If LINALG_FLAVOR has been set then use that library
  find_package(MKL REQUIRED)
  target_include_directories(tango_im.x PUBLIC ${MKL_INCLUDE_DIR})
  target_link_libraries(tango_im.x PUBLIC ${MKL_LIBRARIES})
  target_compile_definitions(tango_im.x PUBLIC USE_MKL)
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
    target_compile_definitions(tango_im.x PUBLIC USE_MKL)
    set(linalg_library_name "MKL")
  else()
    set(linalg_library_name "system default")
  endif()
endif()
message("Using linear algebra library: ${linalg_library_name}")

