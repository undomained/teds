// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "radiative_transfer.h"

#include <common/constants.h>

#include <vector>

auto rtACT(PyObject* /* self */, PyObject* args) -> PyObject*
{
    import_array();
    // Atmospheric gas columns and their cross-sections
    PyArrayObject* np_co2 {};
    PyArrayObject* np_ch4 {};
    PyArrayObject* np_h2o {};
    PyArrayObject* np_xsec_co2 {};
    PyArrayObject* np_xsec_ch4 {};
    PyArrayObject* np_xsec_h2o {};
    // Sentinel 2 albedo B11 band
    PyArrayObject* np_albedo {};
    // Cosines of SZA and VZA
    PyArrayObject* np_mu_sza {};
    PyArrayObject* np_mu_vza {};
    // Solar spectrum
    PyArrayObject* np_sun {};
    // SGM radiances for one scan-line (output)
    PyArrayObject* np_rad {};
    if (!PyArg_ParseTuple(args,
                          "O!O!O!O!O!O!O!O!O!O!O!",
                          &PyArray_Type,
                          &np_co2,
                          &PyArray_Type,
                          &np_ch4,
                          &PyArray_Type,
                          &np_h2o,
                          &PyArray_Type,
                          &np_xsec_co2,
                          &PyArray_Type,
                          &np_xsec_ch4,
                          &PyArray_Type,
                          &np_xsec_h2o,
                          &PyArray_Type,
                          &np_albedo,
                          &PyArray_Type,
                          &np_mu_sza,
                          &PyArray_Type,
                          &np_mu_vza,
                          &PyArray_Type,
                          &np_sun,
                          &PyArray_Type,
                          &np_rad)) {
        return 0;
    }
    double* co2 { reinterpret_cast<double*>(PyArray_DATA(np_co2)) };
    double* ch4 { reinterpret_cast<double*>(PyArray_DATA(np_ch4)) };
    double* h2o { reinterpret_cast<double*>(PyArray_DATA(np_h2o)) };
    double* xsec_co2 { reinterpret_cast<double*>(PyArray_DATA(np_xsec_co2)) };
    double* xsec_ch4 { reinterpret_cast<double*>(PyArray_DATA(np_xsec_ch4)) };
    double* xsec_h2o { reinterpret_cast<double*>(PyArray_DATA(np_xsec_h2o)) };
    double* albedo { reinterpret_cast<double*>(PyArray_DATA(np_albedo)) };
    double* mu_sza { reinterpret_cast<double*>(PyArray_DATA(np_mu_sza)) };
    double* mu_vza { reinterpret_cast<double*>(PyArray_DATA(np_mu_vza)) };
    double* sun { reinterpret_cast<double*>(PyArray_DATA(np_sun)) };
    double* rad { reinterpret_cast<double*>(PyArray_DATA(np_rad)) };
    auto dims { PyArray_DIMS(np_co2) };
    const auto n_act { static_cast<size_t>(dims[0]) };
    const auto n_lay { static_cast<size_t>(dims[1]) };
    dims = PyArray_DIMS(np_xsec_co2);
    const auto n_lbl { static_cast<size_t>(dims[0]) };
    std::vector<double> tau_tot(n_lbl, 0.0);
#pragma omp parallel for firstprivate(tau_tot)
    for (size_t i_act = 0; i_act < n_act; ++i_act) {
        // Python set_opt_depth_species function
        for (size_t i_lbl {}; i_lbl < n_lbl; ++i_lbl) {
            for (size_t i_lay {}; i_lay < n_lay; ++i_lay) {
                const size_t idx_xs { i_lbl * n_lay + i_lay };
                const size_t idx_col { i_act * n_lay + i_lay };
                tau_tot[i_lbl] += (xsec_co2[idx_xs] * co2[idx_col]
                                   + xsec_ch4[idx_xs] * ch4[idx_col]
                                   + xsec_h2o[idx_xs] * h2o[idx_col]);
            }
        }
        for (double& val : tau_tot) {
            val *= 1e-4;
        }
        // Python transmission function
        const double mu_eff { std::abs(1.0 / mu_sza[i_act])
                              + std::abs(1.0 / mu_vza[i_act]) };
        const double fact { mu_sza[i_act] / std::numbers::pi };
        for (size_t i_lbl {}; i_lbl < n_lbl; ++i_lbl) {
            const double exp_tot { std::exp(-tau_tot[i_lbl] * mu_eff) };
            const size_t idx { i_act * n_lbl + i_lbl };
            rad[idx] = sun[i_lbl] * fact * albedo[i_act] * exp_tot;
        }
    }
    return PyLong_FromLong(0);
}
