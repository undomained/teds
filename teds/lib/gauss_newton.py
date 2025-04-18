# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Gauss-Newton optimization."""
import numpy as np
import numpy.typing as npt

from .convolution import Kernel
from .radiative_transfer import OpticAbsProp
from .radiative_transfer import nonscat_fwd_model
from .surface import Surface
from teds.l1l2.types import L2
from teds.sgm.atmosphere import Atmosphere


def gauss_newton(retrieval_init: dict,
                 atm: Atmosphere,
                 optics: OpticAbsProp,
                 wave_lbl: npt.NDArray[np.float64],
                 sun_lbl: npt.NDArray[np.float64],
                 ymeas: npt.NDArray[np.floating],
                 Smeas_diag: npt.NDArray[np.floating],
                 mu0: float,
                 muv: float,
                 isrf: Kernel,
                 timings: dict,
                 l2: L2,
                 i_alt: int,
                 i_act: int) -> None:
    """Nonlinear least squares fit using Gauss-Newton.

    Parameters
    ----------
    retrieval_init
        Initialization dictionary of the least squares fit
    atm
        Atmosphere properties, including trace gas concentrations
    optics
        Optics properties for computing the optical depths
    wave_lbl
        Line-by-line wavelength grid
    sun_lbl
        Solar irradiance on the line-by-line grid
    ymeas
        Radiance spectrum of a given ALT/ACT point
    Smeas_diag
        Radiance spectrum variances
    mu0
        cos(SZA)
    muv
        cos(VZA)
    isrf
        ISRF
    timings
        Collections of measured times. Helps to identify potential
        bottlenecks.
    l2
        Level 2 product (output)
    i_alt
        L2 ALT bin populated with data on this run
    i_act
        L2 ACT bin populated with data on this run

    """
    # Gases to be retrieved
    gas_names = list(retrieval_init['trace_gases'].keys())

    # State vector, to be populated with gas concentrations and albedo
    # coefficients.
    state_vector = np.empty(
        len(gas_names)+len(retrieval_init['albedo_coefficients']))

    # Start by adding gas scaling initial guesses to the state vector
    for i, gas in enumerate(gas_names):
        xcol = sum(atm.get_gas(gas).concentration) / sum(atm.air)
        # Prior scaling parameter
        state_vector[i] = retrieval_init['trace_gases'][gas]['init'] / xcol

    # Next add albedo coefficient initial guesses. The first entries
    # are always gases. Thus, slicing with [len(gas_names):] allows to
    # access the albedo coefficients.
    state_vector[len(gas_names):] = retrieval_init['albedo_coefficients']

    chi_sqrt = []
    n_dof = len(ymeas) - len(state_vector)
    surface = Surface(wave_lbl)
    for iterations in range(retrieval_init['max_iter']):
        for i_gas, gas in enumerate(gas_names):
            atm.get_gas(gas).concentration[:] = (
                state_vector[i_gas]
                * retrieval_init['trace_gases'][gas]['ref_profile'])

        # Calculate surface data
        surface.get_albedo_poly(list(state_vector[len(gas_names):]))

        # Nonscattering forward model
        fwd_rad, derivatives, derivatives_layers = nonscat_fwd_model(
            gas_names,
            len(retrieval_init['albedo_coefficients']),
            isrf,
            sun_lbl,
            atm,
            optics,
            surface,
            mu0,
            muv,
            timings)

        ytilde = ymeas - fwd_rad

        # Calculate Jacobian from convolved quantities
        jacobian = np.empty([len(fwd_rad), len(state_vector)])
        for i, derivative in enumerate(derivatives):
            jacobian[:, i] = derivative

        # Calculate least square solution
        Syinv = np.eye(fwd_rad.size) / Smeas_diag
        # Covariance matrix of the estimated least square solution
        JTJ = np.matmul(jacobian.T, np.matmul(Syinv, jacobian))
        Sx = np.linalg.solve(JTJ, np.eye(JTJ.shape[0]))
        gain = np.matmul(Sx, np.matmul(jacobian.T, Syinv))
        # Least squares solution
        state_delta = np.matmul(gain, ytilde)

        x_lst_precision = np.empty(len(state_delta))
        for m in range(len(state_delta)):
            x_lst_precision[m] = np.sqrt(Sx[m, m])

        state_vector += state_delta

        chi_sqrt.append(((ymeas - fwd_rad)**2 / Smeas_diag).sum() / n_dof)

        if iterations > 2 and (chi_sqrt[-2] - chi_sqrt[-1]
                               < retrieval_init['chi2_lim']):
            l2.converged[i_alt, i_act] = True
            break
        iterations += 1

    l2.chi2[i_alt, i_act] = chi_sqrt[-1]
    l2.number_iter[i_alt, i_act] = iterations

    # Define output product, first update all parameters
    for i_gas, gas in enumerate(gas_names):
        atm.get_gas(gas).concentration[:] = (
            state_vector[i_gas]
            * retrieval_init['trace_gases'][gas]['ref_profile'])

    # Calculate column mixing ratio, precision, and albedo values
    for i_gas, gas in enumerate(gas_names):
        if gas == 'CH4':
            scaling = 1e-9
        else:
            scaling = 1e-6
        mixing_ratio = (
            sum(atm.get_gas(gas).concentration)/sum(atm.air))
        l2.mixing_ratios[gas][i_alt, i_act] = mixing_ratio
        l2.precisions[gas][i_alt, i_act] = (
            mixing_ratio * x_lst_precision[i_gas])
        ref_mixing_ratio = (
            state_vector[i_gas]
            * np.sum(retrieval_init['trace_gases'][gas]['ref_profile'])
            / sum(atm.air) / scaling)
        l2.gains[gas][i_alt, i_act, :] = gain[i_gas] * ref_mixing_ratio
    for i in range(len(retrieval_init['albedo_coefficients'])):
        l2.albedo_coeffs[i, i_alt, i_act] = state_vector[len(gas_names)+i]

    # Calculate column averaging kernels
    n_lay = atm.zlay.size
    for i_gas, gas in enumerate(gas_names):
        col = np.sum(retrieval_init['trace_gases'][gas]['ref_profile'])
        for k_lay in range(n_lay):
            delta_prof = (
                col
                / retrieval_init['trace_gases'][gas]['ref_profile'][k_lay])
            l2.col_avg_kernels[gas][i_alt, i_act, k_lay] = (
                np.sum(gain[i_gas, :] * derivatives_layers[i_gas][:, k_lay])
                * delta_prof)
