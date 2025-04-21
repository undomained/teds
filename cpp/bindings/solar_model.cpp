// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "solar_model.h"

#include <common/solar_model.h>

auto solarModel(PyObject* /* self */, PyObject* args) -> PyObject*
{
    import_array();
    uint32_t tai_seconds {};
    double tai_second_fraction {};
    PyArrayObject* np_q_j2000_ecef {};
    if (!PyArg_ParseTuple(args,
                          "IdO!",
                          &tai_seconds,
                          &tai_second_fraction,
                          &PyArray_Type,
                          &np_q_j2000_ecef)) {
        return 0;
    }
    std::array<double, tango::dims::vec> sun {};
    tango::Quaternion q_j2000_ecef {};
    tango::solarModel(tai_seconds, tai_second_fraction, sun, q_j2000_ecef);
    double* np_q_j2000_ecef_data { reinterpret_cast<double*>(
      PyArray_DATA(np_q_j2000_ecef)) };
    for (int i {}; i < tango::dims::vec; ++i) {
        np_q_j2000_ecef_data[i] = q_j2000_ecef[i];
    }
    np_q_j2000_ecef_data[tango::dims::quat - 1] = q_j2000_ecef.real();
    return PyLong_FromLong(0);
}
