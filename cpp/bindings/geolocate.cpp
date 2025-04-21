// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "geolocate.h"

#include "convert.h"

#include <common/constants.h>
#include <common/geometry.h>

auto geolocate(PyObject* /* self */, PyObject* args) -> PyObject*
{
    import_array();
    char* dem_filename {};
    PyArrayObject* np_los {};
    PyArrayObject* np_tai_seconds {};
    PyArrayObject* np_tai_subsec {};
    PyArrayObject* np_orb_pos {};
    PyArrayObject* np_att_quat {};
    PyArrayObject* np_latitude {};
    PyArrayObject* np_longitude {};
    PyArrayObject* np_height {};
    PyArrayObject* np_vza {};
    PyArrayObject* np_vaa {};
    PyArrayObject* np_sza {};
    PyArrayObject* np_saa {};
    if (!PyArg_ParseTuple(args,
                          "sO!O!O!O!O!O!O!O!O!O!O!O!",
                          &dem_filename,
                          &PyArray_Type,
                          &np_los,
                          &PyArray_Type,
                          &np_tai_seconds,
                          &PyArray_Type,
                          &np_tai_subsec,
                          &PyArray_Type,
                          &np_orb_pos,
                          &PyArray_Type,
                          &np_att_quat,
                          &PyArray_Type,
                          &np_latitude,
                          &PyArray_Type,
                          &np_longitude,
                          &PyArray_Type,
                          &np_height,
                          &PyArray_Type,
                          &np_vza,
                          &PyArray_Type,
                          &np_vaa,
                          &PyArray_Type,
                          &np_sza,
                          &PyArray_Type,
                          &np_saa)) {
        return 0;
    }
    const int n_alt { static_cast<int>(*PyArray_DIMS(np_tai_seconds)) };
    const int n_act { static_cast<int>(*PyArray_DIMS(np_los)) };
    std::vector<double> los(n_act * tango::dims::vec);
    convert2cpp(np_los, los);
    std::vector<uint32_t> tai_seconds(n_alt);
    convert2cpp(np_tai_seconds, tai_seconds);
    std::vector<double> tai_subsec(n_alt);
    convert2cpp(np_tai_subsec, tai_subsec);
    std::vector<double> orb_pos(n_alt * tango::dims::vec);
    convert2cpp(np_orb_pos, orb_pos);
    std::vector<tango::Quaternion> att_quat(n_alt);
    convert2cpp(np_att_quat, att_quat);
    tango::Geometry geo {};
    try {
        tango::geolocate(
          dem_filename, los, tai_seconds, tai_subsec, orb_pos, att_quat, geo);
    } catch (const std::runtime_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    }
    convert2np(geo.lat, np_latitude);
    convert2np(geo.lon, np_longitude);
    convert2np(geo.height, np_height);
    convert2np(geo.vza, np_vza);
    convert2np(geo.vaa, np_vaa);
    convert2np(geo.sza, np_sza);
    convert2np(geo.saa, np_saa);
    return PyLong_FromLong(0);
}
