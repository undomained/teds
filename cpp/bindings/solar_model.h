// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "init.h"

// Using the solar model, return the J2000-to-ECEF rotation quaternion
// using the TAI seconds as input.
auto solarModel(PyObject* /* self */, PyObject* args) -> PyObject*;
