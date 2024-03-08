#==============================================================================
# using the local config to generate the im_config.cfg file
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys
import yaml

def im_configuration(local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code
    
    #modify it accordingly to the yaml file   
    #ckd

    print("Running im_configuration to create config file for IM executable")
    lines = []
    lines.append('# Configuration file for running the instrument model\n\n')

    #==========================main============================================
    lines.append('[main]\n') 
    # Location of the CKD file
    lines.append('ckd_file_in = ' +local_config['ckd_input'] + '\n')
    # Location of SGM output (input spectra)
    lines.append('l1x_input = ' + local_config['rad_input'] + '\n')
    # Location of binning table file
    lines.append('binningtable_filename = ' +local_config['binning_table'] + '\n')
    # Location of L1A product (output)
#    lines.append('l1a_outputfile = ' + local_config['lla_output'] + '\n')
    lines.append('l1a_outputfile = ' + local_config['l1a_file'] + '\n')
    # Which NetCDF group to use from the binning table file.
    # Here it is the same as the binning factor in ACT dimension.
    lines.append('binning_table_id = '+str(local_config['settings_IM']['bin_id']) + '\n')
    # log file path
    lines.append('log_file_path = '+str(local_config['IM_log_path']) + '\n')
    # Exposure time in s
    lines.append('exposure_time = '+str(local_config['settings_IM']['exp_time']) + '\n')
    # co-adding 
    lines.append('nr_coadditions ='+ str(local_config['settings_IM']['co_adding']) + '\n')
    if(local_config['select_images']):
        # The first ALT location to be included in processing (default 0)
        lines.append('image_start = '+str(local_config['first_image'])+'\n')
        # The last ALT location to be included in processing (default all of
        # them). For experimental purposes you may use a smaller number to
        # save time.
        lines.append('# image_end = '+str(local_config['last_image']) + '\n')    
    # enable sub-module of the instrument model (1 = switch on, 0 = switch off)
    # dark current
    lines.append('dark_apply = '+str(local_config['settings_IM']['sw_dark'])+'\n')
    # non-linearity
    lines.append('nonlin_apply = '+str(local_config['settings_IM']['sw_nonlin'])+'\n')
    # spectrometer stray light
    lines.append('stray_apply = '+str(local_config['settings_IM']['sw_stray'])+'\n')
    # pixel response non-uniformity
    lines.append('prnu_apply = '+str(local_config['settings_IM']['sw_prnu'])+'\n')
    # radiometric calibration
    lines.append('rad_apply = '+str(local_config['settings_IM']['sw_rad'])+'\n')
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
    new_config = open(local_config['IM_cfg_path']+'im_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()
    print("DONE running im_configuration to create config file for IM executable")
    return


if __name__ == '__main__' and __package__ is None:

    config = yaml.safe_load(open(sys.argv[1]))
    im_configuration(config)

