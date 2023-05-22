# using the local config to generate the im_config.cfg file

def im_configuration(paths, global_config, local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code

    file_ext     = global_config['profile']+'_'+global_config['run_id']+'.nc\n'
    
    # read in an existing IM config file    
    im_config = open(paths.project+paths.IM_module+'im_config_baseline.cfg')
    lines = im_config.readlines()
    im_config.close()

    #modify it accordingly to the yaml file   
    #ckd
    lines = []
    lines.append('[main]\n') 
    # Location of the CKD file
    lines.append('ckd_file_in = ' + paths.project+paths.data_interface+paths.interface_ckd+local_config['filename'] ['ckd_input']+'.nc\n')
    # Location of SGM output (input spectra)
    lines.append('l1x_input = ' + paths.project+paths.data_interface+paths.interface_sgm + local_config['filename']['sgm_input']+ file_ext)
    # Location of binning table file
    lines.append('binningtable_filename = ' + paths.project+paths.data_interface+paths.interface_ckd + local_config['filename']['binning_table']+'.nc\n')
    # Location of L1A product (output)
    lines.append('l1a_outputfile = ' + paths.project+paths.data_interface+paths.interface_l0 + local_config['filename']['level1a_input']+file_ext)
    # Which NetCDF group to use from the binning table file.
    # Here it is the same as the binning factor in ACT dimension.
    lines.append('binning_table_id = '+str(local_config['settings']['bin_id']) + '\n')
    # Exposure time in s
    lines.append('exposure_time = '+str(local_config['settings']['exp_time']) + '\n')
    # Number of coadditions applied to detector images. A higher number
    # effectively enhances the dynamical range of the count values stored
    # in an integer format.
    lines.append('nr_coadditions ='+str(local_config['settings']['co_adding']) + '\n')
    if(local_config['select_images']):
        # The first ALT location to be included in processing (default 0)
        lines.append('image_start = '+str(local_config['first_image'])+'\n')
        # The last ALT location to be included in processing (default all of
        # them). For experimental purposes you may use a smaller number to
        # save time.
        lines.append('# image_end = '+str(local_config['last_image']) + '\n')    
    #isrf group
    lines.append('[isrf]\n')
    # Whether to convolve input spectra (execute = 1) with the ISRF. If
    # not (execute = 0), the LBL are simply linearly interpolated onto the
    # wavelength grids provided by the CKD.
    lines.append('execute = ' +str(local_config['isrf']['convolution'])+'\n')
    # If doing the ISRF convolution, use this FWHM, in nm, for the Gaussian kernel
    lines.append('fwhm_gauss = ' + str(local_config['isrf']['fwhm']) + '\n')
    #noise model
    lines.append('[noise]\n')
    # Whether to apply noise to the raw detector images (1) or not (0)
    lines.append('noise_apply = '+str(local_config['noise']['switch'])+'\n')
    #Seed for the random noise generator.
    lines.append('seed = '+str(local_config['noise']['seed']) + '\n')
    #level 1b specific group
    lines.append('[l1b]\n')
    
    # write IM config file 
    new_config = open(paths.project+paths.IM_module+'im_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()
    return

