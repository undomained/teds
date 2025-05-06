# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

add_library(bindings SHARED
  bindings.cpp
  convert.cpp
  drivers.cpp
  solar_model.cpp
  geolocate.cpp
  radiative_transfer.cpp)
target_link_libraries(bindings PRIVATE tango_im tango_l1b tango_l2 common)
set_target_properties(bindings PROPERTIES PREFIX "" SUFFIX ".abi3.so")
