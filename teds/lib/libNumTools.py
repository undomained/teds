# collection of numerical tools
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy
import math
from numba import njit

@njit(cache=True)
def convolution(spectrum, isrf, istart, iend):
    #spectral convolution of the spectrum with the isrf
    sh = isrf.shape
    spectrum_conv = np.zeros(sh[0])
    for iwav in range(sh[0]):
        spectrum_conv[iwav] = np.dot(isrf[iwav, istart[iwav]:iend[iwav]], spectrum[istart[iwav]:iend[iwav]])
    return spectrum_conv


def get_isrf_gaussian(parameter, wave_target, wave_input):
    #get a  Gaussian isrf specified by parameter dictionary
    fwhm = parameter['fwhm']
    isrf = {}
    nwave_target = wave_target.size
    nwave_input = wave_input.size
    isrf["isrf"] = np.zeros((nwave_target, nwave_input))
    const = fwhm**2/(4*np.log(2))
    istart = []
    iend = []
    for ii, wmeas in enumerate(wave_target):
        wdiff = wave_input - wmeas
        istart.append(np.argmin(np.abs(wdiff + 1.5*fwhm)))
        iend.append(np.argmin(np.abs(wdiff - 1.5*fwhm)))
        isrf["isrf"][ii, istart[ii]:iend[ii]] = np.exp(-wdiff[istart[ii]:iend[ii]]**2/const)
        isrf["isrf"][ii, :] = isrf["isrf"][ii, :] / np.sum(isrf["isrf"][ii, :])
    isrf["istart"] = np.array(istart)
    isrf["iend"] = np.array(iend)
    return isrf


def get_isrf_generalized_normal(parameter, wave_target, wave_input):
    # get generalized_normal distribution function as isrf
    isrf = {}
    nwave_target = wave_target.size
    nwave_input = wave_input.size
    isrf["isrf"] = np.zeros((nwave_target, nwave_input))
    fwhm = parameter['fwhm']
    bcoeff = parameter['bcoeff']
    const = np.log(2)**bcoeff/(fwhm*np.math.gamma(1+bcoeff))
    istart = []
    iend = []
    for ii, wmeas in enumerate(wave_target):
        wdiff = wave_input - wmeas
        istart.append(np.argmin(np.abs(wdiff + 1.5*parameter['fwhm'])))
        iend.append(np.argmin(np.abs(wdiff - 1.5*parameter['fwhm'])))
        isrf["isrf"][ii, istart[ii]:iend[ii]] = const*2**(-(2*np.abs(wdiff[istart[ii]:iend[ii]])/fwhm)**(1/bcoeff))
        isrf["isrf"][ii, :] = isrf["isrf"][ii, :] / np.sum(isrf["isrf"][ii, :])
    isrf["istart"] = np.array(istart)
    isrf["iend"] = np.array(iend)
    return isrf


def get_isrf(wave_target, wave_input, parameter):
    """Compute kernel and return a convolution function

    Parameters
    ----------
    wave_target : Array (Float)
        Target waves used for convolution
    wave_input : Array(Float)
        wave input
    parameter : Dict
        Dictionary containing settings for convolution

    Return
    --------
    Convolution function which takes a signal as input

    """
    if parameter['type'] == 'Gaussian':
        isrf = get_isrf_gaussian(parameter, wave_target, wave_input)
    elif (parameter['type'] == 'generalized_normal'):
        isrf = get_isrf_generalized_normal(parameter, wave_target, wave_input)

    # convolution function with numba cannot be directly added here
    def conv(spectrum, isrf=isrf):
        return convolution(spectrum, isrf['isrf'], isrf['istart'], isrf['iend'])
    return conv



class isrfct:
    def __init__(self, wave_target, wave_input):

        self.wave_target = wave_target
        self.wave_input = wave_input
        self.isrf = {}

    def get_isrf(self, parameter):

        nwave_target = self.wave_target.size
        nwave_input = self.wave_input.size
        self.isrf['isrf'] = np.zeros((nwave_target, nwave_input))

        if(parameter['type'] == 'Gaussian'):
            const = parameter['fwhm']**2/(4*np.log(2))
            istart = []
            iend = []
            for l, wmeas in enumerate(self.wave_target):
                wdiff = self.wave_input - wmeas
                istart.append(np.argmin(np.abs(wdiff + 1.5*parameter['fwhm'])))
                iend.append(np.argmin(np.abs(wdiff - 1.5*parameter['fwhm'])))

                self.isrf['isrf'][l, istart[l]:iend[l]] = np.exp(-wdiff[istart[l]:iend[l]]**2/const)
                self.isrf['isrf'][l, :] = self.isrf['isrf'][l, :] / np.sum(self.isrf['isrf'][l, :])

            self.isrf['istart'] = istart
            self.isrf['iend'] = iend

        if (parameter['type'] =='generalized_normal'):
            fwhm = parameter['fwhm']
            bcoeff = parameter['bcoeff']
            const = np.log(2)**bcoeff/(fwhm*math.gamma(1+bcoeff))

            istart = []
            iend = []
            for l, wmeas in enumerate(self.wave_target):
                wdiff = self.wave_input - wmeas
                istart.append(np.argmin(np.abs(wdiff + 1.5*parameter['fwhm'])))
                iend.append(np.argmin(np.abs(wdiff - 1.5*parameter['fwhm'])))

                self.isrf['isrf'][l, istart[l]:iend[l]] = const*2**(-(2*np.abs(wdiff[istart[l]:iend[l]])/fwhm)**(1/bcoeff))
                self.isrf['isrf'][l, :] = self.isrf['isrf'][l, :] / np.sum(self.isrf['isrf'][l, :])

            self.isrf['istart'] = istart
            self.isrf['iend'] = iend

        return

    def isrf_convolution(self, spectrum):

        nwave_target = self.wave_target.size
        spectrum_conv = np.empty(nwave_target)

        for iwav in range(nwave_target):
            istart = self.isrf['istart'][iwav]
            iend = self.isrf['iend'][iwav]
            spectrum_conv[iwav] = self.isrf['isrf'][iwav, istart:iend].dot(spectrum[istart:iend])
        return spectrum_conv


def Gaussian2D(size, fwhm_x, fwhm_y, center=None):
    # Make a 2 dimensional gaussian kernel as a product of two Gaussian
    # fwhm are given in units of pixel samples

    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2) / fwhm_x**2) \
        * np.exp(-4*np.log(2) * ((y-y0)**2) / fwhm_y**2)


def getconvolutionparams(kernel_settings, dx, dy):
    """Convolution parameters.

    Parameters
    ----------
    kernel_settings : Dict
        Parameters needed for convolution.
    dx : Float
        Spacing in x.
    dy : Float
        Spacing in y.
    """
    conv_settings = {}
    if (kernel_settings['type'] == '2D Gaussian'):
        fwhm_x = kernel_settings['fwhm_x']
        fwhm_y = kernel_settings['fwhm_y']
        fsize = kernel_settings['size_factor']
        conv_settings['type'] = kernel_settings['type']
        conv_settings['1D kernel extension'] = np.int0(fsize*np.max([fwhm_x, fwhm_y])/np.min([dx, dy]))
        # convert all kernel parameter in units of sampling distance
        conv_settings['fwhm x'] = np.int0(float(fwhm_x)/dx)
        conv_settings['fwhm y'] = np.int0(float(fwhm_y)/dy)
    return conv_settings


def convolution_2d(data, settings):
    # convolve the data array with a kernel defined in settings
    if (settings['type'] == '2D Gaussian'):
        kernel = Gaussian2D(settings['1D kernel extension'], settings['fwhm x'], settings['fwhm y'])

    kernel = kernel / kernel.sum()
    data_conv = scipy.signal.convolve(data, kernel, mode='same')
    return data_conv


def print_attributes(class_object):
    attributes = [attr for attr in dir(class_object)
                  if not attr.startswith('__')]
    print(attributes)
    return

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx], idx

#=========================================================================

class TransformCoords:
    """A class to transform coordinates"""

    def __init__(self, origin):
        """Initialize class based on the Origin 

        Parameters
        ----------
        origin : [lat, lon] Float64
            Origin in lat-lon coordinates
        """
        self.rd = np.pi / 180.0
        phi0 = origin[0] * self.rd
        self.ld0 = origin[1] * self.rd
        self.s_p0 = np.sin(phi0)
        self.c_p0 = np.cos(phi0)
        self.fact = 6371.0
        self.factmts = self.fact * 1000

    def latlon2xykm(self, lat, lon):
        """latlon2xykm Convert from lat lon to xy in km

        Parameters
        ----------
        lat : Matrix/Vector/Float64
            Latitude
        lon : Matrix/Vector/Float64
            Longitude

        Returns
        -------
        x  : Matrix/Vector/Float64
            x in km
        y  : Matrix/Vector/Float64
            y in km
        """
        ld = lon * self.rd
        phi = lat * self.rd
        s_p = np.sin(phi)
        c_p = np.cos(phi)
        ll = ld - self.ld0
        c_l = np.cos(ll)
        s_l = np.sin(ll)
        c_pl = c_p * c_l
        w = np.sqrt(2.0 / (np.maximum(1.0 + self.s_p0 * s_p + self.c_p0 * c_pl, 1.0e-10)))
        x = c_p * s_l * w
        y = (self.c_p0 * s_p - self.s_p0 * c_pl) * w
        return x * self.fact, y * self.fact

    def latlon2xymts(self, lat, lon):
        """latlon2xymts Convert from lat lon to xy in km

        Parameters
        ----------
        lat : Matrix/Vector/Float64
            Latitude
        lon : Matrix/Vector/Float64
            Longitude

        Returns
        -------
        x  : Matrix/Vector/Float64
            x in mts
        y  : Matrix/Vector/Float64
            y in mts
        """
        ld = lon * self.rd
        phi = lat * self.rd
        s_p = np.sin(phi)
        c_p = np.cos(phi)
        ll = ld - self.ld0
        c_l = np.cos(ll)
        s_l = np.sin(ll)
        c_pl = c_p * c_l
        w = np.sqrt(2.0 / (np.maximum(1.0 + self.s_p0 * s_p + self.c_p0 * c_pl, 1.0e-10)))
        x = c_p * s_l * w
        y = (self.c_p0 * s_p - self.s_p0 * c_pl) * w
        return x * self.factmts, y * self.factmts

    def xykm2latlon(self, x1, y1):
        """xykm2latlon Convert from x, y to lat-lon

        Parameters
        ----------
        x1 : Matrix/Vector/Float64
            x coordinate in km
        y1 : Matrix/Vector/Float64
            y coordinate in km

        Returns
        -------
        lat  : Matrix/Vector/Float64
            latitude
        lon  : Matrix/Vector/Float64
            Longitude
        """
        x, y = x1 / self.fact, y1 / self.fact
        p = np.maximum(np.sqrt(x**2 + y**2), 1.0e-10)
        c = 2.0 * np.arcsin(p / 2.0)
        s_c = np.sin(c)
        c_c = np.cos(c)
        phi = np.arcsin(c_c * self.s_p0 + y * s_c * self.c_p0 / p)
        ld = self.ld0 + np.arctan2(x * s_c, (p * self.c_p0 * c_c - y * self.s_p0 * s_c))
        lat = phi / self.rd
        lon = ld / self.rd
        if isinstance(lat, np.ndarray):
            lat[lat > 90.0] -= 180.0
            lat[lat < -90.0] += 180.0
            lon[lon > 180.0] -= 360.0
            lon[lon < -180.0] += 360.0
        else:
            if abs(lat) > 90.0:
                if lat > 0:
                    lat = lat - 180.0
                else:
                    lat = lat + 180.0
            if abs(lon) > 180.0:
                if lon > 0:
                    lon = lon - 360.0
                else:
                    lon = lon + 360.0
        return lat, lon

    def xymts2latlon(self, x1, y1):
        """xymts2latlon Convert from x, y to lat-lon

        Parameters
        ----------
        x1 : Matrix/Vector/Float64
            x coordinate in m
        y1 : Matrix/Vector/Float64
            y coordinate in m

        Returns
        -------
        lat  : Matrix/Vector/Float64
            latitude
        lon  : Matrix/Vector/Float64
            Longitude
        """
        x, y = x1 / self.factmts, y1 / self.factmts
        p = np.maximum(np.sqrt(x**2 + y**2), 1.0e-10)
        c = 2.0 * np.arcsin(p / 2.0)
        s_c = np.sin(c)
        c_c = np.cos(c)
        phi = np.arcsin(c_c * self.s_p0 + y * s_c * self.c_p0 / p)
        ld = self.ld0 + np.arctan2(x * s_c, (p * self.c_p0 * c_c - y * self.s_p0 * s_c))
        lat = phi / self.rd
        lon = ld / self.rd
        if isinstance(lat, np.ndarray):
            lat[lat > 90.0] -= 180.0
            lat[lat < -90.0] += 180.0
            lon[lon > 180.0] -= 360.0
            lon[lon < -180.0] += 360.0
        else:
            if abs(lat) > 90.0:
                if lat > 0:
                    lat = lat - 180.0
                else:
                    lat = lat + 180.0
            if abs(lon) > 180.0:
                if lon > 0:
                    lon = lon - 360.0
                else:
                    lon = lon + 360.0
        return lat, lon
