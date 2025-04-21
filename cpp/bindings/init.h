// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Each C++ extension should be initialized with these definitions

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdexcept>
