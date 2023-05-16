
# main script for Tango's end-to-end simulator

import context
from modules.GM.gm import geometry_module
from modules.SGM.sgm import scene_generation_module
from modules.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from modules.L1L2.l1bl2 import level1b_to_level2_processor

import yaml
import shutil
import sys
import subprocess

def convert_im_config(paths, global_config, local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code

    file_ext     = global_config['profile']+'_'+global_config['run_id']+'.nc\n'
    
    # read in an existing IM config file    
    im_config = open(paths.project+paths.IM_module+'im_config_baseline.cfg')
    lines = im_config.readlines()
    im_config.close()

    #modify it accordingly to the yaml file   
    #ckd
    lines[0] = '[main]\n' 
    #ckd general
    lines[1] = 'ckd_file_in = ' + paths.project+paths.data_interface+paths.interface_ckd+local_config['filename'] ['ckd_input']+'.nc\n'
    #sgm data
    lines[2] = 'l1x_input = ' + paths.project+paths.data_interface+paths.interface_sgm + local_config['filename']['sgm_input']+ file_ext
    #level 0 output
    lines[3] = 'l1a_outputfile = ' + paths.project+paths.data_interface+paths.interface_l0 + local_config['filename']['level0_input']+file_ext
    #ckd binning table
    lines[4] = 'binningtable_filename = ' + paths.project+paths.data_interface+paths.interface_ckd + local_config['filename']['binning_table']+'.nc\n'
    #level0 output
    
    # write IM config file 
    new_config = open(paths.project+paths.IM_module+'im_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()
    return

def convert_l0l1b_config(parhs, global_config, local_config):     
    # using the local and global configs from the yaml files we setup the IM config file of the c++ code
    
    file_ext     = global_config['profile']+'_'+global_config['run_id']+'.nc\n'
    
    # read in an existing l0l1b config file    
    l0l1b_config = open(paths.project+paths.L0L1B_module+'l0l1b_baseline.cfg')
    lines = l0l1b_config.readlines()
    l0l1b_config.close()

    #modify it accordingly to the yaml file   
    #ckd
    lines[0] = '[main]\n' 
    #ckd general
    lines[1] = 'ckd_file_in = ' + paths.project + paths.data_interface+paths.interface_ckd  + local_config['filename'] ['ckd_input']+'.nc\n'
    #ckd binning table
    lines[3] = 'binningtable_filename =' + paths.project + paths.data_interface+paths.interface_ckd  + local_config['filename']['binning_table']+'.nc\n'
    #level 0 output
    lines[7] = 'l1a_files = '+ paths.project + paths.data_interface+paths.interface_l0  + local_config['filename']['level0_input']+file_ext
    #level1b data
    lines[8] = 'outputfile = ' + paths.project + paths.data_interface+paths.interface_l1b  + local_config['filename']['level1b_output']+ file_ext
    #gm input
    lines[9] = 'geometry_file = ' + paths.project + paths.data_interface+paths.interface_gm  + local_config['filename']['gm_input']+ file_ext
    
    # write IM config file 
    new_config = open(paths.project+paths.L0L1B_module+'l0l1b_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()

    return

if __name__ == "__main__":
    # end to end global config file
    
    import paths
    import constants

    global_config = yaml.safe_load(open("./e2es_config.yaml"))

    # ======= geometry module ======================================
    # choose baseline GM config
    if(global_config['setup']['flag_gm']):

        shutil.copyfile(paths.project+ paths.GM_modules+'gm_config_baseline.yaml', \
                        paths.project+ paths.GM_modules+'gm_config.yaml',)

        gm_config = yaml.safe_load(open(paths.project+paths.GM_modules+'gm_config.yaml'))
        geometry_module(paths, global_config, gm_config)

    # ======= scene generator module ===============================
    # choose baseline SGM config
    if(global_config['setup']['flag_sgm']):
        shutil.copyfile(paths.project+paths.SGM_module+'sgm_config_baseline.yaml', \
                        paths.project+paths.SGM_module+'sgm_config.yaml',)
        sgm_config = yaml.safe_load(open(paths.project+paths.SGM_module+ "sgm_config.yaml"))
        scene_generation_module(paths, global_config, sgm_config)

    # ======= The instrument model =================================        
        
    if(global_config['setup']['flag_im']):
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        convert_im_config(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    # ======= The L0L1 pocessor ====================================
    if(global_config['setup']['flag_l0l1b']):
        l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
        convert_l0l1b_config(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

    # ======= The simplified IM and L1B model ======================
    if(global_config['setup']['flag_siml1b']):
        # choose baseline simplified IM and L1B config
        shutil.copyfile(paths.project+ paths.SIML1B_module + 'siml1b_config_baseline.yaml',
                        paths.project+ paths.SIML1B_module +'siml1b_config.yaml',)
        siml1b_config = yaml.safe_load(open(paths.project + paths.SIML1B_module + "siml1b_config.yaml"))
        simplified_instrument_model_and_l1b_processor(paths, global_config, siml1b_config)

    # ======= L1 to L2 processor ===================================
    if(global_config['setup']['flag_l1bl2']):
        # choose baseline L1BL2 config
        shutil.copyfile(paths.project+paths.L1L2_module + 'l1bl2_config_baseline.yaml',
                        paths.project+paths.L1L2_module + 'l1bl2_config.yaml',)
        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config.yaml"))
        level1b_to_level2_processor(paths, global_config, l1bl2_config)
