// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "init.h"

// Geolocation function. Using navigation data as input, generate the
// viewing and solar geometries by modifying the arguments in-place.
auto geolocate(PyObject* /* self */, PyObject* args) -> PyObject*;
