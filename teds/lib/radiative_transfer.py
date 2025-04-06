# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Radiative transfer methods"""
from contextlib import redirect_stdout
from dataclasses import dataclass
from netCDF4 import Dataset
from tqdm import tqdm
from typing import Any
import io
import netCDF4 as nc
import numpy as np
import numpy.typing as npt
import os
import sys
import time

from .hapi import db_begin, fetch_by_ids, absorptionCoefficient_Voigt
from .surface import Surface
from teds import log
from teds.lib.convolution import Kernel
from teds.sgm.atmosphere import Atmosphere
import teds.lib.constants as const


def read_sun_spectrum_TSIS1HSRS(filename: str) -> tuple[
        npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Read sun spectrum TSIS-1 HSRS.

    Downloaded from https://lasp.colorado.edu/lisird/data/tsis1_hsrs,
    Coddington et al., GRL, 2021, https://doi.org/10.1029/2020GL091709

    Parameters
    ----------
    filename
        Path to solar spectrum

    Returns
    -------
    out: dictionary with wavelength [nm], irradiance [W m-2 nm-1]

    """
    ds = nc.Dataset(filename)
    wavelengths = ds['Vacuum Wavelength'][:].data
    return (wavelengths,
            ds['SSI'][:] / (const.hplanck * const.clight) * wavelengths * 1e-9)


class MolecularData:
    """Methods for calculating the absorption cross sections of
    molecular absorbers

    """
    def __init__(self,
                 wave: npt.NDArray[np.float64],
                 xsdbpath: str,
                 hp_ids: list[tuple[str, int]]) -> None:
        """Download line parameters from HITRAN web ressource via the
         hapi tools.

        Populates paths to HITRAN parameter files

        Parameters
        ----------
        wave
            Wavelengths [nm]
        xsdbpath
            Path to location where to store the absorption data
        hp_ids
            list of isotopologue ids, format [(name1, id1),(name2,
            id2) ...]  (see hp.gethelp(hp.ISO_ID))

        """
        self.wave = wave
        # Cross section data
        self.xsdb: dict[Any, Any] = {}
        # Check whether input is in range
        if len(hp_ids) == 0:
            log.error('MolecularData.get_data_HITRAN: provide at least '
                      'one species.')
            sys.exit()

        with redirect_stdout(io.StringIO()):
            db_begin(xsdbpath)

        wv_start = self.wave[0]
        wv_stop = self.wave[-1]
        for id in hp_ids:
            key = '%2.2d' % id[1]
            self.xsdb[key] = {}
            self.xsdb[key]['species'] = id[0]
            # Write 1 file per isotopologue
            self.xsdb[key]['name'] = (
                'ID%2.2d_WV%5.5d-%5.5d' % (id[1], wv_start, wv_stop))
            # Check if data files are already inplace, if not: download
            if (not os.path.exists(
                    os.path.join(xsdbpath, self.xsdb[key]['name']+'.data'))
                and not os.path.exists(
                    os.path.join(xsdbpath, self.xsdb[key]['name']+'.header'))):
                # Wavelength input is [nm], hapi requires wavenumbers [1/cm]
                fetch_by_ids(self.xsdb[key]['name'],
                             [id[1]],
                             1e7 / wv_stop,
                             1e7 / wv_start)


@dataclass
class SpeciesProperties:
    """Absorption properties pertaining to a given species."""
    # Species name, e.g. CO2
    name: str
    # Absorption cross-section
    xsec: npt.NDArray[np.float64]
    # Optical thickness (xsec * rho where rho is concentration at some
    # altitude). Depends on wavelength and altitude.
    tau_alt: npt.NDArray[np.float64]


class OpticAbsProp:
    """Methods to calculate optical scattering and absorption
    properties of single-scattering atmosphere.

    """

    def __init__(self, wave:
                 npt.NDArray[np.float64],
                 zlay: npt.NDArray[np.float64]) -> None:
        # Optical properties. One entry per species.
        self.props: list[SpeciesProperties] = []
        # Total optical thickness
        self.taua = np.empty(())
        # Wavelengths [nm]
        self.wave = wave
        # Vertical height layers, midpoints [m]
        self.zlay = zlay

    def get_prop(self, species: str) -> SpeciesProperties:
        """Return x-section and optical thickness from species name."""
        return list(filter(lambda x: x.name == species, self.props))[0]

    def calc_molec_xsec(
            self, molec_data: MolecularData, atm: Atmosphere) -> None:
        """Calculate molecular absorption cross sections.

        Populates
        ---------
        prop['molec_XX']
            dictionary with optical properties with XXXX HITRAN
            identifier code
        prop['molec_XX']['xsec']
            absorption cross sections [wavelength, nlay] [cm2]

        Parameters
        ----------
        molec_data
            molec_data object
        atm_data
            atmosphere_data object

        """
        nlay = self.zlay.size
        nwave = self.wave.size
        nu_samp = 0.005  # Wavenumber sampling [1/cm] of cross sections
        # Loop over all isotopologues, id = HITRAN global isotopologue ID
        for id in molec_data.xsdb.keys():
            species = molec_data.xsdb[id]['species']
            log.info(f'Calculating absorption cross-section of {species}')
            # Write absorption optical depth [nwave,nlay] in
            # dictionary / per isotopologue
            xsec = np.zeros((nwave, nlay))
            # Check whether absorber type is in the atmospheric data structure

            # Loop over all atmospheric layers
            for ki in tqdm(range(len(atm.zlay))):
                pressure = atm.play[ki]
                temp = atm.tlay[ki]
                # Calculate absorption cross section for layer
                nu, xs = absorptionCoefficient_Voigt(
                    SourceTables=molec_data.xsdb[id]['name'],
                    Environment={'p': pressure / const.PSTD, 'T': temp},
                    WavenumberStep=nu_samp)
                dim_nu = nu.size
                nu_ext = np.insert(nu, 0, nu[0]-nu_samp)
                nu_ext = np.append(nu_ext, nu[dim_nu-1]+nu_samp)
                xs_ext = np.insert(xs, 0, 0.)
                xs_ext = np.append(xs_ext, 0.)
                # Interpolate on wavelength grid provided on input
                xsec[:, ki] = np.interp(
                    self.wave, np.flip(1e7 / nu_ext), np.flip(xs_ext))
            self.props.append(SpeciesProperties(species, xsec, np.empty(())))

    def set_opt_depth_species(
            self, atm: Atmosphere, species_names: list[str]) -> None:
        """Calculate absorption optical depth from the various species
        and combines those as specified.

        Parameters
        ----------
        atm
            atmospheric input and species to be combined to total optical depth

        """
        nlay = self.zlay.size
        nwave = self.wave.size

        # Cross sections are given in cm^2, atmospheric densities in
        # m^2
        conv_f = 1e-4

        for name in species_names:
            prop = self.get_prop(name)
            prop.tau_alt = prop.xsec * atm.get_gas(name).concentration * conv_f

        self.taua = np.zeros((nwave, nlay))
        for name in species_names:
            prop = self.get_prop(name)
            self.taua = self.taua + prop.tau_alt

    def xsec_to_file(self, filename: str) -> None:
        """Save cross-section of each trace gas to NetCDF file."""
        nc = Dataset(filename, 'w')
        nc.title = 'Absorption cross-sections of L2 gases'
        nc.description = 'Each group refers to one trace gas'
        nc.createDimension('layer', len(self.zlay))
        nc.createDimension('wavelength', len(self.wave))
        for prop in self.props:
            grp = nc.createGroup(prop.name)
            var = grp.createVariable('xsec', 'f8', ('wavelength', 'layer'))
            var.long_name = 'absorption cross-section'
            var.units = 'm2'
            var.min_value = 0.0
            var.max_value = 1.0
            var[:] = prop.xsec
        nc.close()

    def xsec_from_file(self, filename: str) -> None:
        """Read cross-sections from NetCDF file."""
        nc = Dataset(filename)
        for species in nc.groups:
            self.props.append(SpeciesProperties(
                species, nc[species]['xsec'][:].data, np.empty(())))


def transmission(sun_lbl: npt.NDArray[np.float64],
                 optics: OpticAbsProp,
                 surface: Surface,
                 mu0: float,
                 muv: float,
                 deriv: bool = False) -> (npt.NDArray[np.float64]
                                          | tuple[npt.NDArray[np.float64],
                                                  npt.NDArray[np.float64],
                                                  npt.NDArray[np.float64]]):
    """Calculate transmission solution given geometry (mu0,muv) using
    matrix algebra

    Parameters
    ----------
    sun_lbl
        Solar irradiance spectrum
    optics
        optic_prop object
    surface
        surface_prop object
    mu0
        cosine of the solar zenith angle
    muv
        cosine of the viewing zenith angle

    Returns
    -------
        Single scattering relative radiance [wavelength] [1/sr]

    """
    if not (0 <= mu0 <= 1 and -1 <= muv <= 1):
        log.error('transmission: input out of range')
        sys.exit(1)

    # Total vertical optical thickness per layer (Delta tau_k) [n_wave,nlay]
    tau_k = optics.taua
    # total optical thickness per spectral bin [n_wave]
    tau_tot = np.sum(tau_k, axis=1)
    mueff = abs(1/mu0) + abs(1/muv)
    fact = mu0 / np.pi
    exp_tot = np.exp(-tau_tot * mueff)

    dev_alb = fact * exp_tot * sun_lbl
    rad_trans = surface.alb * dev_alb

    if deriv:
        # This is the derivative with respect to tau_tot and tau_k
        # because d_tau_tot / d_tau_k = 1.
        dev_tau = -mueff * rad_trans
        return rad_trans, dev_tau, dev_alb
    else:
        return rad_trans


def nonscat_fwd_model(gas_names: list[str],
                      n_albedo: int,
                      isrf: Kernel,
                      sun_lbl: npt.NDArray[np.float64],
                      atm: Atmosphere,
                      optics: OpticAbsProp,
                      surface: Surface,
                      mu0: float,
                      muv: float,
                      timings: dict) -> tuple[npt.NDArray[np.float64],
                                              list[npt.NDArray[np.float64]],
                                              list[npt.NDArray[np.float64]]]:
    """Non-scattering forward model.

    Parameters
    ----------
    gas_names
        Derivatives with respect to these gases will be computed
    n_albedo
        How many albedo coefficients to compute
    isrf
        ISRF
    sun_lbl
        Line-by-line solar irradiance spectrum
    atm
        Atmosphere with trace gas concentrations
    optics
        Optics properties for computing the optical depths
    surface
        Surface properties for recomputing albedo
    mu0
        cos(SZA)
    muv
        cos(VZA)
    timings
        Collections of measured times. Helps to identify potential
        bottlenecks.

    Returns
    -------
        Forward model radiance, derivatives of gas scalings and albedo
        coefficients, and derivative of gas scaling factors per layer
        optical depth.

    """
    time_opt = time.time()
    optics.set_opt_depth_species(atm, gas_names)
    timings['opt'] += time.time() - time_opt

    time_rtm = time.time()
    rad_lbl, dev_tau_lbl, dev_alb_lbl = transmission(
        sun_lbl, optics, surface, mu0, muv, True)
    timings['rtm'] += time.time() - time_rtm

    time_conv = time.time()
    rad = isrf.convolve(rad_lbl)
    derivatives_albedo = []
    derivatives_gases = []
    derivatives_gases_layers = []
    for i in range(n_albedo):
        derivatives_albedo.append(isrf.convolve(dev_alb_lbl * surface.spec**i))
    timings['conv'] += time.time() - time_conv

    time_kern = time.time()
    n_wave = rad.size
    nlay = optics.get_prop(gas_names[0]).tau_alt[0, :].size
    for gas in gas_names:
        # Derivative with respect to a scaling of the total optical depth
        derivatives_gases.append(isrf.convolve(
            np.sum(optics.get_prop(gas).tau_alt, axis=1) * dev_tau_lbl))
        # Derivative with respect to a scaling of the layer optical depth
        layer = np.empty((n_wave, nlay))
        for k_lay in range(optics.get_prop(gas).tau_alt[0, :].size):
            layer[:, k_lay] = isrf.convolve(
                optics.get_prop(gas).tau_alt[:, k_lay] * dev_tau_lbl)
        derivatives_gases_layers.append(layer)
    timings['kern'] += time.time() - time_kern

    return (rad,
            derivatives_gases+derivatives_albedo,
            derivatives_gases_layers)
