// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions to convert between Python and C++ arrays

#pragma once

#include "init.h"

#include <vector>

namespace tango {

class Quaternion;

} // namespace tango

// Convert a Numpy array to a C++ vector
template <typename T>
auto convert2cpp(const PyArrayObject* np_data, std::vector<T>& data) -> void
{
    T* np_ptr { reinterpret_cast<T*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data[i] = np_ptr[i];
    }
}

// Convert a Numpy array containing quaternions to a C++ vector
auto convert2cpp(const PyArrayObject* np_data,
                 std::vector<tango::Quaternion>& data) -> void;

// Convert a C++ vector to a Numpy array
auto convert2np(const std::vector<double>& data,
                PyArrayObject* np_data) -> void;
