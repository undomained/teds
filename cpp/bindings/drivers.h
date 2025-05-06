// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Define Python-C++ bindings for all modules written in C++

#pragma once

#include "init.h"

auto runIM(PyObject* /* self */, PyObject* args) -> PyObject*;
auto runL1B(PyObject* /* self */, PyObject* args) -> PyObject*;
auto runL2(PyObject* /* self */, PyObject* args) -> PyObject*;
