# using the yaml file create the corresponding l0l1b_config.cfg to run the L0L1B processor
def l0l1b_configuration(paths, global_config, local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code
    
    file_ext     = global_config['profile']+'_'+global_config['run_id']+'.nc\n'
    
    # read in an existing l0l1b config file    
    l0l1b_config = open(paths.project+paths.L0L1B_module+'l0l1b_baseline.cfg')
    lines = l0l1b_config.readlines()
    l0l1b_config.close()

    #modify it accordingly to the yaml file   
    #ckd
    lines = []
    #man group
    lines.append('[main]\n') 
    # Process name. In the present implementation this is always l1b.
    lines.append('process = l1b\n')
    #location of ckd file
    lines.append('ckd_file_in = ' + paths.project + paths.data_interface+paths.interface_ckd  + local_config['filename'] ['ckd_input']+'.nc\n')
    #location of ckd binning table
    lines.append('binningtable_filename =' + paths.project + paths.data_interface+paths.interface_ckd  + local_config['filename']['binning_table']+'.nc\n')
    #l1b group
    lines.append('[l1b]\n') 
    # Iteration limit for the stray light deconvolution procedure. Set to
    # 0 turn off the stray light correction.
    lines.append('stray_van_cittert_steps = '+str(local_config['settings']['van_cittert_steps'])+'\n')
    lines.append('stray_interpolating = '+str(local_config['settings']['stray_interpolating'])+'\n')
    # Whether to determine geolocation. Use 0 for now.
    lines.append('geolocation = '+str(local_config['settings']['geolocation']) +'\n')
    #Location of the nL1A product (input detector images)
    lines.append('l1a_files = '+ paths.project + paths.data_interface+paths.interface_l0  + str(local_config['filename']['level0_input'])+file_ext)
    #Location of the level1b data
    lines.append('outputfile = ' + paths.project + paths.data_interface+paths.interface_l1b  + local_config['filename']['level1b_output']+ file_ext)
    #location of gm input
    lines.append('geometry_file = ' + paths.project + paths.data_interface+paths.interface_gm  + local_config['filename']['gm_input']+ file_ext)
    
    # write IM config file 
    new_config = open(paths.project+paths.L0L1B_module+'l0l1b_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()

    return
