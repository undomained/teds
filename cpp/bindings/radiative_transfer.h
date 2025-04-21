// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "init.h"

// Radiative transfer for one scan-line (across all ACT points)
auto rtACT(PyObject* /* self */, PyObject* args) -> PyObject*;
