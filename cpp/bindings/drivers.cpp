// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "drivers.h"

#include "convert.h"

#include <im/driver_im.h>
#include <im/settings_im.h>
#include <iostream>
#include <l1al1b/driver_l1b.h>
#include <l1al1b/settings_l1b.h>
#include <l1l2/driver_l2.h>
#include <l1l2/settings_l2.h>

// Read Python argument into a dictionary, parse that to a string and
// initialize the settings class instance. Return true if
// successful. Otherwise, if the Python argument was None, print all
// configuration parameters for the given module and exit. This
// function is applicable to all C++ modules of TEDS.
static auto settingsInit(PyObject* args, tango::Settings& settings) -> bool
{
    PyObject* conf_py {};
    if (!PyArg_ParseTuple(args, "|O", &conf_py)) {
        return false;
    }
    std::string conf {};
    if (conf_py != nullptr && conf_py != Py_None) {
        if (PyDict_Check(conf_py)) {
            PyObject* repr { PyObject_Repr(conf_py) };
            conf = PyUnicode_AsUTF8(repr);
            Py_DECREF(repr);
        } else {
            PyErr_SetString(PyExc_TypeError, "Argument must be a dict or None");
            return false;
        }
    }
    if (conf.empty()) {
        std::cout << "Stopping because module was called without an argument.\n"
                     "Valid configuration parameters for this module:\n\n"
                  << settings.c_str() << '\n';
        return false;
    }
    settings.yamlLoad(conf);
    settings.init();
    return true;
}

auto runIM(PyObject* /* self */, PyObject* args) -> PyObject*
{
    tango::SettingsIM settings {};
    if (settingsInit(args, settings)) {
        tango::driverIM(settings);
    }
    return PyLong_FromLong(0);
}

auto runL1B(PyObject* /* self */, PyObject* args) -> PyObject*
{
    tango::SettingsL1B settings {};
    if (settingsInit(args, settings)) {
        tango::driverL1B(settings);
    }
    return PyLong_FromLong(0);
}

auto runL2(PyObject* /* self */, PyObject* args) -> PyObject*
{
    tango::SettingsL2 settings {};
    if (settingsInit(args, settings)) {
        tango::driverL2(settings);
    }
    return PyLong_FromLong(0);
}
