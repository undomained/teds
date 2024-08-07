"""
Instrument Spectral Response Function
"""
import numpy as np
from teds.CKD.ckd_generation.generator_class import *
import os
from netCDF4 import Dataset
import csv
from scipy.interpolate import CubicSpline


def import_data(conf): # Get the data and combine
    print("[isrf] >> Importing data")
    isrf_nc = Dataset(conf['paths']['dir_input']+'isrf_data.nc', 'w', format="NETCDF4")
    data_dir = conf['paths']['dir_external'] + 'TANGO_Nitro_011_ISRF_SonyTolerancedMp2S_LowResSampling/'
    data_files = os.listdir(data_dir)

    wl, ap = [], []  #wavelengths, act_pos
    # Get files, as well as wavelength and actpos data (from filenames...)
    for fi, f in enumerate(data_files):
        parts = f.split('_')
        wl.append(float(parts[-2][:-2])) # wavelength in nm
        ap.append(int(parts[-1][:-10])) # across track position 

    # Combine data in dataset
    wavelengths, actpos = np.unique(wl), np.unique(ap)
    for i, file in enumerate(data_files):
        with open(data_dir + file, 'r', newline='\n') as f:
            rows = csv.reader(f, delimiter = ',')
            filedata = np.array([row for row in rows])
            headers = filedata[0, :]
            filedata = filedata[1:, :]

            if i == 0:  # initialize data array
                data = np.zeros((len(wavelengths), len(actpos), len(filedata), len(filedata[0])))
                isrf_nc.createDimension('wavelength', len(wavelengths))
                var = isrf_nc.createVariable('wavelength', 'float32', 'wavelength')
                var[:] = wavelengths
                var.units = 'nm'

                isrf_nc.createDimension('actpos', len(actpos))
                var = isrf_nc.createVariable('actpos', 'u1', 'actpos')
                var[:] = actpos

                isrf_nc.createDimension('datapoints_isrf', len(filedata[:,0]))
                for h, head in enumerate(headers[:2]):
                    var = isrf_nc.createVariable(headers[h], 'float32', 'datapoints_isrf')
                    units = 'nm' if 'nm' in head else 'um'
                    var.units = units
                    var[:] = filedata[:,h]
            
            p = np.flatnonzero(wavelengths == wl[i])
            q = np.flatnonzero(actpos == ap[i])
            data[p,q,:,:] = filedata

    var = isrf_nc.createVariable(headers[2], 'float32', ['wavelength', 'actpos', 'datapoints_isrf'])
    var[:] = data[:,:,:,2]
    isrf_nc.close()


def calculate(conf): # Interpolate and Calculate ISRF
    filepath = conf['paths']['dir_input']+'isrf_data.nc'
    if not os.path.isfile(filepath):
        import_data(conf)
    isrf_nc = Dataset(filepath, 'a', format="NETCDF4")
    # Now interpolate and add to dataset
    isrf_data = isrf_nc['ISRF_Mp2S']
    actpos = isrf_nc['actpos']
    actpos_lbl = ['actpos0 (nadir?)', 'actpos1 (full swath?)']
    wavelengths = np.array(isrf_nc['wavelength'])
    distance_physical = np.array(isrf_nc['Distance_Physicalspace_um'])
    distance_spectral = np.array(isrf_nc['Distance_Wavelengthspace_nm'])

    spectral_dim_name = 'spectral_detector_pixels' 
    spatial_dim_name  = 'spatial_samples_per_image'
    isrf_samples_dim_name = 'isrf_samples'
    X = conf['dimensions'][spectral_dim_name]  # samples in spectral direction
    Y = conf['dimensions'][spatial_dim_name]  # samples in spatial direction
    N_data = len(distance_physical)            # input number of samples of ISRF shape
    N = conf['RADIOMETRIC'][isrf_samples_dim_name]  # output number of samples of ISRF shape
    ix_picks = np.linspace(start = 0, stop = N_data-1, num = N, dtype = int)

    # init arrays for interpolation in spectral direction
    intpdata_wl = np.zeros((len(actpos), X, N))  
    wl_new = np.linspace(np.min(wavelengths), np.max(wavelengths), X)
    ap_new = np.linspace(np.min(actpos), np.max(actpos), Y)

    # Interpolate ISRFs in spectral space
    for ap in  actpos:
        data = isrf_data[:,ap,:]
        
        for i_int, i_dat in enumerate(ix_picks):
            spl = CubicSpline(wavelengths, data[:,i_dat])
            intpdata_wl[ap,:,i_int] = spl(wl_new)

        for wl in range(len(wavelengths)):
            if wl < len(wavelengths) - 1:
                intwl = (wavelengths[wl] + wavelengths[wl+1])/2
                i = np.flatnonzero(wl_new >= intwl)[0]

    # init array with interpolated data in spectral direction
    isrf_interpolated = np.zeros((Y, X, N), dtype = 'float32')  #float32 to reduce size

    cntr = 0
    print('[isrf] >> Calculating ISRF matrix... may take a moment')
    for i, wl in enumerate(wl_new):
        progress = ((i+1)/X*100)
        if progress >= cntr:
            print('[isrf] >> {:.0f}% - {:.1f}nm'.format(progress, wl))
            cntr += 10
        for j in range(N):
            spl = CubicSpline(actpos, intpdata_wl[:,i,j])
            isrf_interpolated[:,i,j] = spl(ap_new)

    # Write to isrf input file
    isrf_nc.createGroup('interpolated')
    isrf_nc['interpolated'].createDimension(spectral_dim_name, X)
    isrf_nc['interpolated'].createDimension(spatial_dim_name, Y)
    isrf_nc['interpolated'].createDimension(isrf_samples_dim_name, N)
    varspec = isrf_nc['interpolated'].createVariable(spectral_dim_name, 'float32', spectral_dim_name)
    varspat = isrf_nc['interpolated'].createVariable(spatial_dim_name, 'float32', spatial_dim_name)
    varsamp_physical = isrf_nc['interpolated'].createVariable('physical_distance', 'float32', isrf_samples_dim_name)
    varsamp_spectral = isrf_nc['interpolated'].createVariable('spectral_distance', 'float32', isrf_samples_dim_name)
    varisrf = isrf_nc['interpolated'].createVariable('isrf_interpolated', 'float32', [spatial_dim_name, spectral_dim_name, 'isrf_samples'])
    varspec.units = 'nm'
    varspat.units = 'um'
    varsamp_physical.units = 'um'
    varsamp_spectral.units = 'nm'
    varisrf.units = ''
    varspec[:] = wl_new
    varspat[:] = ap_new
    varsamp_physical[:] = distance_physical[ix_picks]
    varsamp_spectral[:] = distance_spectral[ix_picks]
    varisrf[:] = isrf_interpolated
    isrf_nc.close()


def generate(conf):
    gen = ckd_generator()
    gen.dim_names = ['spatial_samples_per_image', 'spectral_detector_pixels', 'isrf_samples']
    gen.dtype = 'float32'

    filepath = conf['paths']['dir_input'] + 'isrf_data.nc'
    if not os.path.isfile(filepath):
        calculate(conf)

    isrf_nc = Dataset(filepath, 'r', format="NETCDF4")
    try:
        isrf = isrf_nc['interpolated/isrf_interpolated']
    except:
        # If data not found, recalculate
        isrf_nc.close()
        calculate(conf)
        isrf_nc = Dataset(filepath, 'r', format="NETCDF4")
        isrf = isrf_nc['interpolated/isrf_interpolated']
        
    gen.data = isrf
    gen.attr_names = ['source', 'comment']
    gen.attr_vals = ['TANGO_Nitro_011_ISRF_SonyTolerancedMp2S_LowResSampling (TNO)',
                     '[!] Make sure the isrf_sample increments agree with the increments in the line-by-line spectra for correction convolution']

    return gen

