// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions to convert between Python and C++ arrays

#pragma once

#include "init.h"

#include <common/constants.h>
#include <common/eigen.h>
#include <vector>

// Convert Numpy array to C++ vector
template <typename T>
auto convert2cpp(const PyArrayObject* np_data, std::vector<T>& data) -> void
{
    T* np_ptr { reinterpret_cast<T*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data[i] = np_ptr[i];
    }
}

// Convert Numpy array to Eigen array
auto convert2cpp(const PyArrayObject* np_data,
                 ArrayXNd<tango::dims::vec>& data) -> void;

// Convert Numpy quaternion array to Eigen quaternion array
auto convert2cpp(const PyArrayObject* np_data,
                 std::vector<Eigen::Quaterniond>& data) -> void;

// Convert an Eigen array to a Numpy array
auto convert2np(const ArrayXXd& data, PyArrayObject* np_data) -> void;
