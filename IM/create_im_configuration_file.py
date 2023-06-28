#==============================================================================
# using the local config to generate the im_config.cfg file
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys
import numpy as np
import netCDF4 as nc

def im_configuration(locations, local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code
    
    #modify it accordingly to the yaml file   
    #ckd
    lines = []
    #==========================main============================================
    lines.append('[main]\n') 
    # Location of the CKD file
    lines.append('ckd_file_in = ' +locations['ckd_input'] + '\n')
    # Location of SGM output (input spectra)
    lines.append('l1x_input = ' + locations['rad_input'] + '\n')
    # Location of binning table file
    lines.append('binningtable_filename = ' +locations['binning_table'] + '\n')
    # Location of L1A product (output)
    lines.append('l1a_outputfile = ' + locations['output'] + '\n')
    # Which NetCDF group to use from the binning table file.
    # Here it is the same as the binning factor in ACT dimension.
    lines.append('binning_table_id = '+str(local_config['settings']['bin_id']) + '\n')
    # Exposure time in s
    lines.append('exposure_time = '+str(local_config['settings']['exp_time']) + '\n')
    #========Determine optimal exposure time (sec) and co-adding factor from LUTs =========================

    gm_data = nc.Dataset(locations['gm_input'], mode='r')
    sza_avg = np.mean(gm_data['sza'][:, :])
    #LUT for optimal exposure and coadding factors    
    fmc_lut = np.arange(1,6)
    data    = np.genfromtxt(locations['ckd_expo_time'])
    sza_lut = data[:,0]
    exposure_time_lut = data[:,1:]
    data    = np.genfromtxt(locations['ckd_coadding'])
    coad_fact_lut     = np.int16(data[:,1:])
    #find right indices
    ind_sza = np.argmin(np.abs(sza_lut-sza_avg))
    ind_fmc = np.argmin(np.abs(fmc_lut-local_config['settings']['fmc']))
    lines.append('exposure_time = '+str(exposure_time_lut[ind_sza,ind_fmc]) + '\n')
    # Number of coadditions applied to detector images. A higher number
    # effectively enhances the dynamical range of the count values stored
    # in an integer format.
    lines.append('nr_coadditions ='+str(coad_fact_lut[ind_sza,ind_fmc]) + '\n')
    if(local_config['select_images']):
        # The first ALT location to be included in processing (default 0)
        lines.append('image_start = '+str(local_config['first_image'])+'\n')
        # The last ALT location to be included in processing (default all of
        # them). For experimental purposes you may use a smaller number to
        # save time.
        lines.append('# image_end = '+str(local_config['last_image']) + '\n')    
    # enable sub-module of the instrument model (1 = switch on, 0 = switch off)
    # dark current
    lines.append('dark_apply = '+str(local_config['settings']['sw_dark'])+'\n')
    # non-linearity
    lines.append('nonlin_apply = '+str(local_config['settings']['sw_nonlin'])+'\n')
    # spectrometer stray light
    lines.append('stray_apply = '+str(local_config['settings']['sw_stray'])+'\n')
    # pixel response non-uniformity
    lines.append('prnu_apply = '+str(local_config['settings']['sw_prnu'])+'\n')
    # radiometric calibration
    lines.append('rad_apply = '+str(local_config['settings']['sw_rad'])+'\n')
    #==========================isrf============================================
    lines.append('[isrf]\n')
    # Whether to convolve input spectra (execute = 1) with the ISRF. If
    # not (execute = 0), the LBL are simply linearly interpolated onto the
    # wavelength grids provided by the CKD.
    lines.append('execute = ' +str(local_config['isrf']['convolution'])+'\n')
    # If doing the ISRF convolution, use this FWHM, in nm, for the Gaussian kernel
    lines.append('fwhm_gauss = ' + str(local_config['isrf']['fwhm']) + '\n')
    #========================noise model=======================================
    lines.append('[noise]\n')
    # Whether to apply noise to the raw detector images (1) or not (0)
    lines.append('noise_apply = '+str(local_config['noise']['switch'])+'\n')
    #Seed for the random noise generator.
    lines.append('seed = '+str(local_config['noise']['seed']) + '\n')
    
    # write IM config file 
    new_config = open(locations['IM_path']+'im_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()
    return

