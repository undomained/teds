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
from teds.l1l2.types import RefProfiles
from teds.sgm.atmosphere import Atmosphere


def gauss_newton(config: dict,
                 ref_profiles: RefProfiles,
                 atm: Atmosphere,
                 optics: OpticAbsProp,
                 wave_lbl: npt.NDArray[np.float64],
                 sun_lbl: npt.NDArray[np.float64],
                 sun: npt.NDArray[np.floating],
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
    gas_names = config['retrieval']['gases']
    n_gas = len(gas_names)
    state_vector_names = {}
    for i, gas_name in enumerate(gas_names):
        state_vector_names[gas_name] = i
    n_state = n_gas

    n_albedo = config['retrieval']['n_albedos']
    for i in range(n_albedo):
        state_vector_names['albedo_' + str(i)] = i
    n_state += n_albedo

    if config['retrieval']['do_shift']:
        state_vector_names['shift'] = n_state
        n_state += 1

    # State vector, to be populated with gas concentrations and albedo
    # coefficients.
    state_vector = np.empty(n_state)

    # Start by adding gas scaling initial guesses to the state vector
    for i, gas in enumerate(gas_names):
        xcol = sum(atm.get_gas(gas).concentration) / sum(atm.air)
        # Prior scaling parameter
        state_vector[i] = ref_profiles.initial[gas] / xcol

    # Next add albedo coefficient initial guesses. The first entries
    # are always gases. Thus, slicing with [n_gas:n_gas+n_albedo]
    # allows to access the albedo coefficients.
    state_vector[n_gas:n_gas+n_albedo] = 0
    # Derive first guess albedo from maximum reflectance
    idx = np.argmax(ymeas)
    state_vector[n_gas] = ymeas[idx] / sun[idx] * np.pi / mu0

    if config['retrieval']['do_shift']:
        state_vector[-1] = 0

    chi2 = []
    n_dof = len(ymeas) - len(state_vector)
    surface = Surface(wave_lbl)
    for iterations in range(config['retrieval']['max_iter']):
        # Update all fitted variables using the current state vector
        for i_gas, gas in enumerate(gas_names):
            atm.get_gas(gas).concentration[:] = (
                state_vector[i_gas] * ref_profiles.gases[gas])
        surface.get_albedo_poly(list(state_vector[n_gas:n_gas+n_albedo]))
        if config['retrieval']['do_shift']:
            isrf.regenerate(state_vector[n_gas+n_albedo])

        # Nonscattering forward model
        fwd_rad, derivatives, derivatives_layers = nonscat_fwd_model(
            gas_names,
            n_albedo,
            config['retrieval']['do_shift'],
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
        JTJ = np.matmul(jacobian.T, np.matmul(Syinv, jacobian))
        Sx = np.linalg.solve(JTJ, np.eye(JTJ.shape[0]))
        gain = np.matmul(Sx, np.matmul(jacobian.T, Syinv))
        state_delta = np.matmul(gain, ytilde)

        precisions = np.empty(len(state_delta))
        for m in range(len(state_delta)):
            precisions[m] = np.sqrt(Sx[m, m])

        state_vector += state_delta

        chi2.append(((ymeas - fwd_rad)**2 / Smeas_diag).sum() / n_dof)

        if iterations > 2 and (chi2[-2] - chi2[-1]
                               < config['retrieval']['chi2_lim']):
            l2.converged[i_alt, i_act] = True
            break

    l2.chi2[i_alt, i_act] = chi2[-1]
    l2.iterations[i_alt, i_act] = iterations + 1

    # Define output product, first update all parameters
    for i_gas, gas in enumerate(gas_names):
        atm.get_gas(gas).concentration[:] = (
            state_vector[i_gas] * ref_profiles.gases[gas])

    # Calculate column mixing ratio, precision, and albedo values
    for i_gas, gas in enumerate(gas_names):
        if gas == 'CH4':
            scaling = 1e-9
        else:
            scaling = 1e-6
        mixing_ratio = (
            sum(atm.get_gas(gas).concentration)/sum(atm.air))
        l2.mixing_ratios[gas][i_alt, i_act] = mixing_ratio
        l2.precisions[gas][i_alt, i_act] = mixing_ratio * precisions[i_gas]
        ref_mixing_ratio = (state_vector[i_gas]
                            * np.sum(ref_profiles.gases[gas])
                            / sum(atm.air) / scaling)
        l2.gains[gas][i_alt, i_act, :] = gain[i_gas] * ref_mixing_ratio
    l2.albedo0[i_alt, i_act] = state_vector[n_gas]

    # Spectral shift
    if config['retrieval']['do_shift']:
        l2.spec_shift[i_alt, i_act] = state_vector[n_gas + n_albedo]
    else:
        l2.spec_shift[i_alt, i_act] = 0

    # Calculate column averaging kernels
    n_lay = atm.zlay.size
    for i_gas, gas in enumerate(gas_names):
        col = np.sum(ref_profiles.gases[gas])
        for k_lay in range(n_lay):
            delta_prof = col / ref_profiles.gases[gas][k_lay]
            l2.col_avg_kernels[gas][i_alt, i_act, k_lay] = (
                np.dot(gain[i_gas, :], derivatives_layers[i_gas][:, k_lay])
                * delta_prof)
