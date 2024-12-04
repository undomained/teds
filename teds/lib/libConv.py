import sys
import numpy as np
import netCDF4 as nc
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import datetime
from tqdm import tqdm
import pandas as pd

from teds.lib.libWrite import writevariablefromname
from teds import log

def convolve_isrf(xls_file, wvl_signal, n_rows, signal):
    '''
    Read ISRF xls file
    Fit supergaussian (flat top)
    Interpolate to provided signal wavelengths and n_rows
    Calculate ISRF shape and convolve with signal
    '''

    # read ISRF excel
    df_sim = pd.read_excel(xls_file, "SimulationLogic", skiprows = 3, nrows = 181)
    df_sim.columns = ['ix', 'wavelength', 'actpos']
    df_sim['wavelength'] *= 1000 # convert wavelength to nm
    df_wl = pd.read_excel(xls_file, "WavelengthDifferenceSampling", skiprows = 4, nrows = 1029, usecols= "A:FY")
    df_isrf = pd.read_excel(xls_file, "ISRF", skiprows = 5, nrows = 1029, usecols= "H:GF")

    sim_indices = df_sim.iloc[:,0]
    wavelengths = np.unique(df_sim.iloc[:,1])
    actpos = np.unique(df_sim.iloc[:,2])
    
    # fit supergaussian(x, x0, w, n)
    popt = np.zeros((len(sim_indices), 3))
    significant_indices = np.zeros((len(sim_indices)))
    for sim_ix in sim_indices:
        wldif = np.array(df_wl[sim_ix])  # nm
        isrf = np.array(df_isrf[sim_ix]) # AU
        
        p0 = [0, 0.6, 10]
        bmin = [-0.3, 0.05, 1]
        bmax = [0.3, 5, 100]
        popt[sim_ix], pcov = curve_fit(supergaussian, wldif, isrf, p0, bounds = (bmin, bmax))
        
        thr = 0.001
        mask = isrf > thr
        significant_indices[sim_ix] = np.max(np.abs(wldif[mask]))
    
    max_wl_diff_to_compute_isrf_at = np.max(significant_indices) # 0.61

    df_sim["x0"] = popt[:, 0]
    df_sim["w"] = popt[:, 1]
    df_sim["n"] = popt[:, 2]

    
    # parameters are smooth so can be interpolated for ACT and wavelength dimensions of input spectra
    n_wvl = len(wvl_signal)
    params_intermediate = np.zeros((len(actpos), n_wvl, 3))
    params = np.zeros((n_rows, n_wvl, 3))
    # interpolate wavelengths
    for i, ap in enumerate(actpos):
        mask = df_sim.iloc[:,2] == ap
        df_filt = df_sim[mask]
        for k, param in enumerate(["x0", "w", "n"]):
            wl = np.array(df_filt["wavelength"])
            p = np.array(df_filt[param])
            spl = CubicSpline(wl, p)
            params_intermediate[i, :, k] = spl(wvl_signal)
    
    # interpolate act
    new_actpos = np.linspace(actpos[0], actpos[-1], n_rows)
    for i, wl in enumerate(wvl_signal):
        mask = df_sim.iloc[:,1] == wl
        for k, param in enumerate(["x0", "w", "n"]):
            spl = CubicSpline(actpos, params_intermediate[:, i, k])
            params[:,i,k] = spl(new_actpos)

    # supergaussian parameters [row, wvl]
    x0 = params[:,:,0]
    w = params[:,:,1]
    n = params[:,:,2]

    # round up to nearest 0.1 nm for wvl diff range
    max_wl_diff = np.true_divide(np.ceil(max_wl_diff_to_compute_isrf_at * 10**1), 10**1)

    # create delta wvl
    mean_inc_input = np.diff(wvl_signal)[0]

    max_ix_wing = (max_wl_diff / mean_inc_input) + 1; 
    dwvl = np.zeros( int((2 * max_ix_wing + 1)))
    for i in range(len(dwvl)):
        dwvl[i] = (i - max_ix_wing) * mean_inc_input

    n_samples = len(dwvl)

    # apply convolution
    signal_conv = np.zeros((n_rows,n_wvl))
    for i_row in range(n_rows):
        for i_wvl in range(n_wvl):

            # left bound of dwl
            i_dwvl_0 = round(max(max_ix_wing - i_wvl, 0))
            i_dwvl_1 = round(min(n_samples - 1 , n_wvl - i_wvl + max_ix_wing - 1))

            # Calculate in-bound part of ISRF
            this_dwvl = dwvl[i_dwvl_0:i_dwvl_1+1]
            this_isrf = supergaussian(this_dwvl, x0[i_row,i_wvl], w[i_row,i_wvl], n[i_row,i_wvl])
            
            # Normalization factor      
            norm_inv = 1/ np.sum(this_isrf)

            # Left and right bounds of input spectrum,
            i_input_0 = round(max(i_wvl - max_ix_wing, 0))
            i_input_1 = round(min(i_wvl + max_ix_wing, n_wvl - 1))
            
            # Calculate in-bound part of input spectrum
            if signal.ndim == 1:
                this_signal = signal[i_input_0:i_input_1+1]
            elif signal.ndim == 2:
                this_signal = signal[i_row,i_input_0:i_input_1+1]

            # Carry out convolution
            for k in range(len(this_isrf)):
                signal_conv[i_row,i_wvl] += norm_inv * this_isrf[k] * this_signal[k]

    return signal_conv

def supergaussian(x, x0, w, n):
    return np.exp(-(np.abs(x-x0)/w)**(2*n))

def conv_rad_gaussian(sgm_rad_file, fwhm):
    # convolve radiance with Gaussian ISRF

    with nc.Dataset(sgm_rad_file) as f:

        rad = f['radiance'][:,:,:] # [alt, act, spectral_bins] - "photons / ( sr nm m2 s)"
        wvl = f['wavelength'][:] # [spectral_bins] - nm

    stepsize = wvl[1]-wvl[0]
    fwhm_step = fwhm/stepsize
    sigma = fwhm_step /np.sqrt(8*np.log(2))
    convolved_rad = gaussian_filter1d(rad, sigma, axis=-1)

    file_out = sgm_rad_file.replace('.nc','_convolved.nc')
    # open file
    with nc.Dataset(file_out, mode='w') as output_conv_rad:
        output_conv_rad.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        output_conv_rad.comment = f'Convolved radiance with Gaussian FWHM {fwhm} nm'
        output_conv_rad.createDimension('along_track', rad.shape[0])     
        output_conv_rad.createDimension('across_track', rad.shape[1])     
        output_conv_rad.createDimension('wavelength', rad.shape[2])

        rad_out = output_conv_rad.createVariable('radiance', float, ('along_track', 'across_track', 'wavelength'))
        wvl_out = output_conv_rad.createVariable('wavelength', float, ('wavelength'))

        rad_out[:,:,:] = convolved_rad
        wvl_out[:] = wvl


    return file_out


def conv_irr(sgm_rad_file, mode, fwhm=None, isrf_file=None):
    '''
    convolve irradiance

    mode = ['Gaussian', 'ISRF']
    fwhm: float: nm
    isrf_file: str: xls file
    '''

    log.info(f'Convolving irradiance from SGM with {mode}')

    with nc.Dataset(sgm_rad_file) as f:

        irr = f['solar_irradiance'][:] # [spectral_bins] - "photons / (nm m2 s)"
        wvl = f['wavelength'][:] # [spectral_bins] - nm
        nrows = f.dimensions['across_track'].size

    stepsize = wvl[1]-wvl[0]

    if mode == 'Gaussian':
        fwhm_step = fwhm/stepsize
        sigma = fwhm_step /np.sqrt(8*np.log(2))
        convolved_irr = gaussian_filter1d(irr, sigma)
        comment = f'Convolved irradiance with Gaussian FWHM {fwhm} nm'
        dims = ('wavelength',)

    elif mode == 'ISRF':
        convolved_irr = convolve_isrf(isrf_file, wvl, nrows, irr)
        comment = f'Convolved irradiance with ISRF'
        dims = ('across_track','wavelength')

    file_out = sgm_rad_file.replace('.nc','_conv_irr.nc')

    with nc.Dataset(file_out, mode='w') as output_conv_irr:
        
        output_conv_irr.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        output_conv_irr.comment = comment 

        output_conv_irr.createDimension('wavelength', len(irr))     # spectral axis

        if 'across_track' in dims:
            output_conv_irr.createDimension('across_track', nrows) 

        # wavelength
        _ = writevariablefromname(output_conv_irr, 'wavelength', ('wavelength',), wvl)
        # solar irradiance
        _ = writevariablefromname(output_conv_irr, 'solarirradiance', dims, convolved_irr)

    return file_out



def conv_rad(sgm_rad_file, mode, fwhm=None, isrf_file=None):
    '''
    convolve radiance

    mode = ['Gaussian', 'ISRF']
    fwhm: float: nm
    isrf_file: str: xls file
    '''

    log.info(f'Convolving radiance from SGM with {mode}')

    with nc.Dataset(sgm_rad_file) as f:

        rad = f['radiance'][:,:,:] # [alt, act, spectral_bins] - "photons / ( sr nm m2 s)"
        wvl = f['wavelength'][:] # [spectral_bins] - nm
        nrows = f.dimensions['across_track'].size
        nalt = f.dimensions['along_track'].size

    stepsize = wvl[1]-wvl[0]

    if mode == 'Gaussian':
        fwhm_step = fwhm/stepsize
        sigma = fwhm_step /np.sqrt(8*np.log(2))
        convolved_rad = gaussian_filter1d(rad, sigma)
        comment = f'Convolved radiance with Gaussian FWHM {fwhm} nm'

    elif mode == 'ISRF':
        convolved_rad = np.zeros_like(rad)
        for ialt in tqdm(range(nalt)):
            convolved_rad[ialt,:,:] = convolve_isrf(isrf_file, wvl, nrows, rad[ialt,:,:])
        
        comment = f'Convolved radiance with ISRF'

    file_out = sgm_rad_file.replace('.nc','_conv_rad.nc')

    with nc.Dataset(file_out, mode='w') as output_conv_rad:
        
        output_conv_rad.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        output_conv_rad.comment = comment 

        output_conv_rad.createDimension('along_track', nalt)     
        output_conv_rad.createDimension('across_track',nrows)     
        output_conv_rad.createDimension('wavelength', len(wvl))

        rad_out = output_conv_rad.createVariable('radiance', float, ('along_track', 'across_track', 'wavelength'))
        wvl_out = output_conv_rad.createVariable('wavelength', float, ('wavelength'))

        rad_out[:,:,:] = convolved_rad
        wvl_out[:] = wvl

    return file_out
