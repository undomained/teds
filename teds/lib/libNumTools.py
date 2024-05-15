# collection of numerical tools
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy
import math
from numba import njit
from scipy.interpolate import RegularGridInterpolator, griddata
from tqdm import tqdm

class Emptyclass:
    pass

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
#=========================================================================

# N-D linear interpolation

def ndim_lin_interpol_get_indices_weights(xp, x):
    
    ndims = len(xp)

    xf = []
    for i in range(ndims):
        xf.append(x[i].flatten())

    nvals = xf[0].shape[0]

    # nearest to grid points and weight in each dimension 
    points = np.zeros((2,ndims,nvals), dtype=int)
    weights = np.zeros((2,ndims,nvals)) * np.nan

    for idim in range(ndims):
        sign = np.sign(xp[idim][1] - xp[idim][0])

        f_idx = np.interp( sign*xf[idim][:], sign*xp[idim][:], np.arange(xp[idim].shape[0]) )

        p = np.floor( f_idx ).astype(int)
        p[p<0] = 0
        p[p>=xp[idim].shape[0]] = xp[idim].shape[0] - 1
        points[0,idim,:] = p[:]

        p = np.ceil( f_idx ).astype(int)
        p[p<0] = 0
        p[p>=xp[idim].shape[0]] = xp[idim].shape[0] - 1
        points[1,idim,:] = p[:]

        w = 1. - (f_idx % 1)
        w[w<0.] = 0.0
        w[w>1.] = 1.0
        weights[0,idim,:] = w
        weights[1,idim,:] = 1. - w

    # combine the points to 2**nmdim indices and corresponding weights
    n = 2**ndims

    p = np.zeros((n, ndims, nvals), dtype=int)
    w = np.ones((n,nvals))

    for idim in range(ndims):
        idx  = np.mod(( np.arange(n) / 2**idim).astype(int),2) 
        p[:,idim,:] = points[:,idim,:][idx,:]
        w = w * weights[:,idim,:][idx,:]

    idx = []
    for i in range( p.shape[1]):
        idx.append( p[:,i,:].flatten().tolist() )
    idx = tuple(idx)

    return idx, w, x[0].shape

def ndim_lin_interpol_get_values(A, idx, w, shape):
    a = np.sum( np.reshape( A[idx], w.shape ) * w[:,:], axis=0) / np.sum( w[:,:], axis=0)
    return np.reshape(a, shape)


def ndim_lin_interpol(xp, A, x):

    ndims = A.ndim
    
    points = np.zeros((ndims,2), dtype=np.int)
    weights = np.zeros((ndims,2))
    
    for idim in range(ndims):
        i0 = np.argmin( np.abs(x[idim] - xp[idim]) )

        i1 = i0 + np.sign(x[idim] - xp[idim][i0]).astype(int)
        
        if i1>xp[idim].shape[0]-1:
            i1=i0-1
        if i1<0:
            i1=i0+1
        
        xdiff = np.float( xp[idim][i1] - xp[idim][i0] )
        if np.abs(xdiff) < 2*np.finfo(float).eps:
            w0 = 1.
        else: 
            w0 =  1. - np.float( x[idim] - xp[idim][i0] ) / xdiff
        
        w1 = 1. - w0
        
        if w1<0:
            w1=0
            w0=1
        if w1>1:
            w1=0
            w0=1 
        
        points[idim,0] = i0
        points[idim,1] = i1
        
        weights[idim,0] = w0
        weights[idim,1] = w1


    n = 2**ndims
    p = np.zeros((n, ndims), dtype=int)
    w = np.ones(n)
    
    for idim in range(ndims):  
        idx  = np.mod(( np.arange(n) / 2**idim).astype(int),2) 
        p[:,idim] = points[idim,:][idx]
        w = w * weights[idim,:][idx]

    idx = tuple(p.T.tolist())
    a = np.sum( A[idx] * w)

    return a, idx, w


def ndim_lin_interpol_evaluate(A, idx, w):
    return np.sum( A[idx] * w)

    
# if __name__ == "__main__":
#     import scipy.interpolate

#     dims = np.array([4, 5, 6, 3, 5, 4, 8])
#     ndim = dims.shape[0]
    
#     A = np.arange(dims.prod()).reshape(dims).astype(np.float)
    
#     xp = []
#     for i in range(ndim):
#         xp_n = np.arange(dims[i])
#         xp.append(xp_n)
        
#     x = np.array([1.1, 2.7, 4.5, 0.1, 3., 0.6, 2.4])
    
#     res = ndim_lin_interpol(xp, A, x)
    
#     res2 = scipy.interpolate.interpn(xp, A, x)
    
#     print(res[0])
    
def interpolate_data_irregular(indata, gm_data, gases):
    """Interpolate irregular data in lat-lon coordinates.

    Parameters
    ----------
    indata : Class 
        Class containing meterological data
    gm_data : Dict
        Dictionary containing the parameters from geometry module.
    config : List
        Contains lits of gases to be processed
    """
    # create a new container for SGM meteo data
    outdata = Emptyclass()
    outdata.__setattr__("lat", gm_data.lat)
    outdata.__setattr__("lon", gm_data.lon)
    outdata.__setattr__("zlev", indata.zlev)
    outdata.__setattr__("zlay", indata.zlay)

    print('Interpolating data to GM mesh...')

    dim_alt, dim_act = outdata.lat.shape   # dimensions
    dxdy = np.column_stack((indata.lat.ravel(), indata.lon.ravel()))

    # Interpolate values to GM grid
    albedo = griddata(dxdy, indata.albedo_conv.ravel(), (outdata.lat, outdata.lon), fill_value=0.0)
    for gas in gases:
        conv_gas = indata.__getattribute__(gas+"_conv")
        interpdata = np.zeros([dim_alt, dim_act, indata.zlay.size])
        for iz in tqdm(range(indata.zlay.size)):
            interpdata[:, :, iz] = griddata(dxdy, conv_gas[:,:, iz].ravel(), (outdata.lat, outdata.lon), fill_value=0.0)
        outdata.__setattr__(gas, interpdata)
    print('                     ...done')
    return albedo, outdata

def interpolate_data_regular(indata, gm_data, gases):
    """Interpolate data on regular cartesian grid.

    Meteo data is in regular x-y grid and this is used to interpolate
    the values to x-y of gm grid. This is implementation is faster
    compare to interpolating data in lat-lon grid.

    Parameters
    ----------
    indata : Class
        Meteo data
    gm_data : Dict
        Parameters from geometry module
    gases : List
        List of gases to be processed

    Returns
    -------
    albedo : Matrix
        Albedo on gm grid
    outdata: Class
        Meteo data on gm grid
    """
    print('Interpolating data to GM mesh...')
    outdata = Emptyclass()
    outdata.__setattr__("lat", gm_data.lat)
    outdata.__setattr__("lon", gm_data.lon)

    dim_alt, dim_act = outdata.lat.shape   # dimensions
    dim_lay = indata.zlay.shape[2]
    dim_lev = indata.zlev.shape[2]
    
    # Interpolate albedo to GM grid
    fa = RegularGridInterpolator((indata.ypos, indata.xpos), indata.albedo_conv, 
                                 bounds_error=False, fill_value=0.0)
    albedo = fa((gm_data.ypos, gm_data.xpos))
    for gas in gases:
        interpdata = np.zeros([dim_alt, dim_act, dim_lay])
        conv_gas = indata.__getattribute__("dcol_"+gas+"_conv")
        for iz in tqdm(range(dim_lay)):
            gasmin = np.min(conv_gas[:, :, iz])
            fa = RegularGridInterpolator((indata.ypos, indata.xpos), conv_gas[:, :, iz], 
                                         bounds_error=False, fill_value=gasmin)
            interpdata[:, :, iz] = fa((gm_data.ypos, gm_data.xpos))
        outdata.__setattr__(gas.upper(), interpdata)

    #Here, we make a shortcut using a equally vertically gridded atmosphere 
    #So we duplicate the vertical grid!
    dum = np.tile(indata.zlay[0,0,:], dim_alt*dim_act)
    outdata.__setattr__("zlay", np.reshape(dum,(dim_alt,dim_act,dim_lay)))

    dum = np.tile(indata.col_air[0,0], dim_alt*dim_act)
    outdata.__setattr__("air", np.reshape(dum,(dim_alt,dim_act)))

    dum = np.tile(indata.zlev[0,0,:], dim_alt*dim_act)
    outdata.__setattr__("zlev", np.reshape(dum,(dim_alt,dim_act,dim_lev)))
     
    print('                             ...done')
    return albedo, outdata

def convolvedata(data, config):

    """Convolve meteo and albedo data.

    Parameters
    ----------
    meteodata : Class
        Meteo data
    config : Dict
        Dict containing configuration parameters.

    """
    print('Convolution data...')
    
    dx = np.mean(np.diff(data.xpos))
    dy = np.mean(np.diff(data.ypos))
    
    conv_settings = getconvolutionparams(config['kernel_parameter'], dx, dy)
    
    # convolution of albedo   
    data.__setattr__("albedo_conv", convolution_2d(data.albedo, conv_settings))
    

    for gas in config['conv_gases']:
        concgas = data.__getattribute__("dcol_"+gas)
        conv_gas = np.zeros_like(concgas)
        for iz in tqdm(range(data.zlay[0,0,:].size)):
            # convolution
            conv_gas[:, :, iz] = convolution_2d(concgas[:, :, iz], conv_settings)
            data.__setattr__("dcol_"+gas+"_conv", conv_gas)

    print('                     ...done')
    
    return data

