# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Collection of numerical tools."""
import numpy as np
import scipy


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
        conv_settings['1D kernel extension'] = np.intp(
            fsize*np.max([fwhm_x, fwhm_y])/np.min([dx, dy]))
        # convert all kernel parameter in units of sampling distance
        conv_settings['fwhm x'] = float(fwhm_x)/dx
        conv_settings['fwhm y'] = float(fwhm_y)/dy
    return conv_settings


def convolution_2d(data, settings):
    # convolve the data array with a kernel defined in settings
    if (settings['type'] == '2D Gaussian'):
        kernel = Gaussian2D(settings['1D kernel extension'],
                            settings['fwhm x'],
                            settings['fwhm y'])

    kernel = kernel / kernel.sum()
    return scipy.signal.convolve(data, kernel, mode='same')


# =========================================================================

# N-D linear interpolation

def ndim_lin_interpol_get_indices_weights(xp, x):

    ndims = len(xp)

    xf = []
    for i in range(ndims):
        xf.append(x[i].flatten())

    nvals = xf[0].shape[0]

    # nearest to grid points and weight in each dimension
    points = np.zeros((2, ndims, nvals), dtype=int)
    weights = np.zeros((2, ndims, nvals)) * np.nan

    for idim in range(ndims):
        sign = np.sign(xp[idim][1] - xp[idim][0])

        f_idx = np.interp(sign*xf[idim][:], sign*xp[idim][:],
                          np.arange(xp[idim].shape[0]))

        p = np.floor(f_idx).astype(int)
        p[p < 0] = 0
        p[p >= xp[idim].shape[0]] = xp[idim].shape[0] - 1
        points[0, idim, :] = p[:]

        p = np.ceil(f_idx).astype(int)
        p[p < 0] = 0
        p[p >= xp[idim].shape[0]] = xp[idim].shape[0] - 1
        points[1, idim, :] = p[:]

        w = 1. - (f_idx % 1)
        w[w < 0.] = 0.0
        w[w > 1.] = 1.0
        weights[0, idim, :] = w
        weights[1, idim, :] = 1. - w

    # combine the points to 2**nmdim indices and corresponding weights
    n = 2**ndims

    p = np.zeros((n, ndims, nvals), dtype=int)
    w = np.ones((n, nvals))

    for idim in range(ndims):
        idx = np.mod((np.arange(n) / 2**idim).astype(int), 2)
        p[:, idim, :] = points[:, idim, :][idx, :]
        w = w * weights[:, idim, :][idx, :]

    idx = []
    for i in range(p.shape[1]):
        idx.append(p[:, i, :].flatten().tolist())
    idx = tuple(idx)

    return idx, w, x[0].shape


def ndim_lin_interpol_get_values(A, idx, w, shape):
    a = (np.sum(np.reshape(A[idx], w.shape) * w[:, :], axis=0)
         / np.sum(w[:, :], axis=0))
    return np.reshape(a, shape)


def ndim_lin_interpol(xp, A, x):

    ndims = A.ndim

    points = np.zeros((ndims, 2), dtype=np.int)
    weights = np.zeros((ndims, 2))

    for idim in range(ndims):
        i0 = np.argmin(np.abs(x[idim] - xp[idim]))

        i1 = i0 + np.sign(x[idim] - xp[idim][i0]).astype(int)

        if i1 > xp[idim].shape[0]-1:
            i1 = i0-1
        if i1 < 0:
            i1 = i0+1

        xdiff = np.float(xp[idim][i1] - xp[idim][i0])
        if np.abs(xdiff) < 2*np.finfo(float).eps:
            w0 = 1.
        else:
            w0 = 1. - np.float(x[idim] - xp[idim][i0]) / xdiff

        w1 = 1. - w0

        if w1 < 0:
            w1 = 0
            w0 = 1
        if w1 > 1:
            w1 = 0
            w0 = 1

        points[idim, 0] = i0
        points[idim, 1] = i1

        weights[idim, 0] = w0
        weights[idim, 1] = w1

    n = 2**ndims
    p = np.zeros((n, ndims), dtype=int)
    w = np.ones(n)

    for idim in range(ndims):
        idx = np.mod((np.arange(n) / 2**idim).astype(int), 2)
        p[:, idim] = points[idim, :][idx]
        w = w * weights[idim, :][idx]

    idx = tuple(p.T.tolist())
    a = np.sum(A[idx] * w)

    return a, idx, w
