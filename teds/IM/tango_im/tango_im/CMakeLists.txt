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
