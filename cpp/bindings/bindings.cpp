// Python interface to TEDS C++ components

#include "drivers.h"
#include "geolocate.h"
#include "radiative_transfer.h"
#include "solar_model.h"

static PyMethodDef bindings_methods[] = {
    { "run_instrument_model",
      runIM,
      METH_VARARGS,
      "Driver function to run the instrument model" },
    { "run_l1al1b",
      runL1B,
      METH_VARARGS,
      "Driver function to run the L1B processor" },
    { "level1b_to_level2_processor",
      runL2,
      METH_VARARGS,
      "Driver function to run the L2 processor" },
    { "solar_model",
      solarModel,
      METH_VARARGS,
      "Return sun vector and ECI-to-ECEF quaternion using the solar model" },
    { "geolocate",
      geolocate,
      METH_VARARGS,
      "Return solar and viewing geometries" },
    { "rt_act",
      rtACT,
      METH_VARARGS,
      "Radiative transfer calculation for one scan-line" },
    { 0, 0, 0, 0 }
};

static struct PyModuleDef bindingsmodule = {
    PyModuleDef_HEAD_INIT,
    "bindings",
    "Python interface to TEDS C++ modules",
    -1,
    bindings_methods
};

PyMODINIT_FUNC PyInit_bindings(void)
{
    auto py_object { PyModule_Create(&bindingsmodule) };
    return py_object;
}
