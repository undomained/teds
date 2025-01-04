# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
from collections import namedtuple
import logging
import numpy as np

from .exceptions import ProcessError

_logger = logging.getLogger(__name__)


RefractiveModCoefficients = namedtuple(
    "RefractiveModCoefficients",
    ['wave_nm', 'mr_BC', 'mi_BC', 'mr_IO', 'mi_IO',
     'mr_DUST', 'mi_DUST', 'mr_OC', 'mi_OC', 'mr_H2O',
     'mi_H2O'])


def replaceD(s):
    s = s.replace(b'D', b'E')
    return s.replace(b'd', b'e')


class RefractiveIndex(object):
    def __init__(self, path):

        self.refractive_path = path
        self.pathBC = self.refractive_path + 'refractive_index_BC.dat'
        self.pathIO = self.refractive_path + 'refractive_index_INORG.dat'
        self.pathDUST = self.refractive_path + 'refractive_index_DUST.dat'
        self.pathOC = self.refractive_path + 'refractive_index_ORG.dat'
        self.pathH2O = self.refractive_path + 'refractive_index_water.dat'

        self.refractive_index_parameters = [
            "BC", "INORG", "DUST", "ORG", 'H2O']
        self.dataBC = None
        self.dataIO = None
        self.dataDUST = None
        self.dataOC = None
        self.dataH2O = None

    def read(self):

        _logger.info("Reading refractive index tables.")

        self.dataBC = np.loadtxt(
            self.pathBC,
            skiprows=1,
            delimiter=None,
            converters={0: replaceD, 1: replaceD, 2: replaceD})
        self.dataIO = np.loadtxt(
            self.pathIO,
            skiprows=1,
            delimiter=None,
            converters={0: replaceD, 1: replaceD, 2: replaceD})
        self.dataDUST = np.loadtxt(
            self.pathDUST,
            skiprows=1,
            delimiter=None,
            converters={0: replaceD, 1: replaceD, 2: replaceD})
        self.dataOC = np.loadtxt(
            self.pathOC,
            skiprows=1,
            delimiter=None,
            converters={0: replaceD, 1: replaceD, 2: replaceD})
        self.dataH2O = np.loadtxt(
            self.pathH2O,
            skiprows=1,
            delimiter=None,
            converters={0: replaceD, 1: replaceD, 2: replaceD})

    def get_coefficients_for_wavelength(self, wavelength_nm):

        if self.dataBC is None:
            self.read()

        try:
            # Lookup the wavelengths.
            wavelengths = np.array(self.dataBC[:, 0])
            # Find the index of the closest wavelength.
            diff_array = np.abs(wavelengths - wavelength_nm)
            indexBC = diff_array.argmin()
            mr_BC = self.dataBC[indexBC, 1].item()
            mi_BC = self.dataBC[indexBC, 2].item()

            wavelengths = np.array(self.dataIO[:, 0])
            diff_array = np.abs(wavelengths - wavelength_nm)
            indexIO = diff_array.argmin()
            mr_IO = self.dataIO[indexIO, 1].item()
            mi_IO = self.dataIO[indexIO, 2].item()

            wavelengths = np.array(self.dataDUST[:, 0])
            diff_array = np.abs(wavelengths - wavelength_nm)
            indexDUST = diff_array.argmin()
            mr_DUST = self.dataDUST[indexDUST, 1].item()
            mi_DUST = self.dataDUST[indexDUST, 2].item()

            wavelengths = np.array(self.dataOC[:, 0])
            diff_array = np.abs(wavelengths - wavelength_nm)
            indexOC = diff_array.argmin()
            mr_OC = self.dataOC[indexOC, 1].item()
            mi_OC = self.dataOC[indexOC, 2].item()

            wavelengths = np.array(self.dataH2O[:, 0])
            diff_array = np.abs(wavelengths - wavelength_nm)
            indexH2O = diff_array.argmin()
            mr_H2O = self.dataH2O[indexH2O, 1].item()
            mi_H2O = self.dataH2O[indexH2O, 2].item()

        except Exception as e:
            raise ProcessError(
                "Error in determining the refractive index coefficients on a "
                "certain wavelength: {}".format(e))

        return RefractiveModCoefficients(wave_nm=wavelength_nm,
                                         mr_BC=mr_BC,
                                         mi_BC=mi_BC,
                                         mr_IO=mr_IO,
                                         mi_IO=mi_IO,
                                         mr_DUST=mr_DUST,
                                         mi_DUST=mi_DUST,
                                         mr_OC=mr_OC,
                                         mi_OC=mi_OC,
                                         mr_H2O=mr_H2O,
                                         mi_H2O=mi_H2O)

    def collocate(self, julday_orbit, lat, lon, npixels):
        pass

    def write_to_group(self, group, start=0, end=None):

        if self.dataBC is None:
            self.read()

        _logger.info("Writing data species")

        if group.name == 'aerosol':

            nspecies = 6
            # Two times same data (also in echam 'parts')
            data_species = [self.dataBC,  self.dataDUST, self.dataH2O,
                            self.dataOC, self.dataIO, self.dataIO]

            # for refractive index

            for i in range(nspecies):
                cols = data_species[i]

                wavelength = np.array(cols[:, 0])
                rri = np.array(cols[:, 1])
                iri = np.array(cols[:, 2])

                var = 'wavelength_species' + str(i + 1)
                group.variables[var][:] = wavelength
                var = 'RRI_species' + str(i + 1)
                group.variables[var][:] = rri
                var = 'IRI_species' + str(i + 1)
                group.variables[var][:] = iri

        if group.name == 'cloud':
            cols = self.dataH2O
            wavelength = np.array(cols[:, 0])
            rri = np.array(cols[:, 1])
            iri = np.array(cols[:, 2])

            var = 'wavelength_water_refr'
            group.variables[var][:] = wavelength
            var = 'RRI_water_refr'
            group.variables[var][:] = rri
            var = 'IRI_water_refr'
            group.variables[var][:] = iri
