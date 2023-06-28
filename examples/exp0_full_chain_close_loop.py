#!/bin/env python
# main script for Tango's end-to-end simulator

# import context
from GM.gm import geometry_module
from SGM.sgm import scene_generation_module
from IM.create_im_configuration_file import im_configuration
from L0L1B.create_l01b_configuration_file import l0l1b_configuration
from SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from L1L2.l1bl2 import level1b_to_level2_processor

import yaml
import shutil
import sys
import subprocess

if __name__ == "__main__":
    # end to end global config file
    
    import paths
    import constants

    global_config = yaml.safe_load(open("./exp0_config.yaml"))

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
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    # ======= The L0L1 pocessor ====================================
    if(global_config['setup']['flag_l0l1b']):
        l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
        l0l1b_configuration(paths, global_config, l0l1b_config)     
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
