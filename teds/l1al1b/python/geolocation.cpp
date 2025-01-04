// Python interface to L1B geolocation routines.

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdexcept>

#include "../tango_l1b/geometry.h"
#include "../tango_l1b/solar_model.h"

// Using the solar model, return the J2000-to-ECEF rotation quaternion
// using the TAI seconds as input.
static PyObject* solarModel(PyObject* /* self */, PyObject* args)
{
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

// Convert a Numpy array to a C++ vector
template <typename T>
static auto convert2cpp(const PyArrayObject* np_data, std::vector<T>& data)
{
    T* np_ptr { reinterpret_cast<T*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data[i] = np_ptr[i];
    }
}

// Convert a Numpy array containing quaternions to a C++ vector
static auto convert2cpp(const PyArrayObject* np_data,
                        std::vector<tango::Quaternion>& data)
{
    double* np_ptr { reinterpret_cast<double*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data[i] = { np_ptr[i * tango::dims::quat + 0],
                    np_ptr[i * tango::dims::quat + 1],
                    np_ptr[i * tango::dims::quat + 2],
                    np_ptr[i * tango::dims::quat + 3] };
    }
}

// Convert a C++ vector to a Numpy array
static auto convert2np(const std::vector<double>& data, PyArrayObject* np_data)
{
    double* np_ptr { reinterpret_cast<double*>(PyArray_DATA(np_data)) };
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        np_ptr[i] = data[i];
    }
}

// Geolocation function. Using navigation data as input, generate the
// viewing and solar geometries by modifying the arguments in-place.
static PyObject* geolocate(PyObject* /* self */, PyObject* args)
{
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

static PyMethodDef geolocation_methods[] = {
    { "solar_model",
      solarModel,
      METH_VARARGS,
      "Return sun vector and ECI-to-ECEF quaternion using the solar model" },
    { "geolocate",
      geolocate,
      METH_VARARGS,
      "Return solar and viewing geometries" },
    { 0, 0, 0, 0 }
};

static struct PyModuleDef geomodule = {
    PyModuleDef_HEAD_INIT,
    "geolocation",
    "Python interface to the L1B processor geolocation routines",
    -1,
    geolocation_methods
};

PyMODINIT_FUNC PyInit_geolocation(void)
{
    auto py_object { PyModule_Create(&geomodule) };
    import_array();
    return py_object;
}
