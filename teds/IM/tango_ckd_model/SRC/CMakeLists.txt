target_sources(ckdmodel PRIVATE
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

get_target_property(_sources ckdmodel SOURCES)
add_custom_target(gitfile
  COMMAND ${CMAKE_COMMAND}
  -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
  -P ${PROJECT_SOURCE_DIR}/cmake/modules/get_git_info.cmake)
add_dependencies(ckdmodel gitfile)

target_include_directories(ckdmodel PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR})

# Link against a linear algebra library
if (LINALG_FLAVOR MATCHES "mkl")
  # If LINALG_FLAVOR has been set then use that library
  find_package(MKL REQUIRED)
  target_include_directories(ckdmodel PUBLIC ${MKL_INCLUDE_DIR})
  target_link_libraries(ckdmodel PUBLIC ${MKL_LIBRARIES})
  target_compile_definitions(ckdmodel PUBLIC USE_MKL)
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
    target_compile_definitions(ckdmodel PUBLIC USE_MKL)
    set(linalg_library_name "MKL")
  else()
    set(linalg_library_name "system default")
  endif()
endif()
message("Using linear algebra library: ${linalg_library_name}")

