// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "convert.h"

#include <common/quaternion.h>

auto convert2cpp(const PyArrayObject* np_data,
                 std::vector<tango::Quaternion>& data) -> void
{
    double* np_ptr { reinterpret_cast<double*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data[i] = { np_ptr[i * tango::dims::quat + 0],
                    np_ptr[i * tango::dims::quat + 1],
                    np_ptr[i * tango::dims::quat + 2],
                    np_ptr[i * tango::dims::quat + 3] };
    }
}

auto convert2np(const std::vector<double>& data, PyArrayObject* np_data) -> void
{
    double* np_ptr { reinterpret_cast<double*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        np_ptr[i] = data[i];
    }
}
