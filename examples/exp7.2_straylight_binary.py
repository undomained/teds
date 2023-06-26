
# exp7.2: stray light analysis for a binary transition scene and with spectral dependence

import context
from end_to_end.GM.gm import geometry_module
from end_to_end.SGM.sgm import scene_generation_module
from end_to_end.IM.create_im_configuration_file import im_configuration
from end_to_end.L0L1B.create_l01b_configuration_file import l0l1b_configuration
from end_to_end.L0L1B.pixel_mask import determine_pixel_mask
from end_to_end.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from end_to_end.L1L2.l1bl2 import level1b_to_level2_processor

import yaml
import shutil
import sys
import subprocess

if __name__ == "__main__":
    # end to end global config file
    
    import paths
    import constants
    
    global_config = yaml.safe_load(open("./exp7.2_straylight_binary.yaml"))

    # ======= geometry module ======================================
    # choose baseline GM config
    if(global_config['setup']['flag_gm']):

        shutil.copyfile(paths.project+ paths.GM_module+'gm_config_baseline.yaml', \
                        paths.project+ paths.GM_module+'gm_config.yaml',)

        gm_config = yaml.safe_load(open(paths.project+paths.GM_module+'gm_config.yaml'))
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
        
        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_exp7.2_21kernel_b1_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_exp7.2_1kernel_b1_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 0
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_exp7.2_wo_strayl_b1_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)
    # ======= The L0L1 pocessor ====================================
    if(global_config['setup']['flag_l0l1b']):

        l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
        l0l1b_config['settings']['van_cittert_steps']   = 4
        l0l1b_config['filename']['level0'] ='Tango_Carbon_l0_exp7.2_21kernel_b1_'                   #level 1a input
        l0l1b_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_21kernel_corrected_b1_'
        l0l1b_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
        l0l1b_configuration(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

        l0l1b_config['settings']['van_cittert_steps'] = 0                 #switch off stray light
        l0l1b_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_21kernel_not_corrected_b1_'
        l0l1b_configuration(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

        l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
        l0l1b_config['settings']['van_cittert_steps'] = 4
        l0l1b_config['filename']['level0']  = 'Tango_Carbon_l0_exp7.2_1kernel_b1_'                   #level 1a input
        l0l1b_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_1kernel_corrected_b1_'
        l0l1b_config['filename']['ckd']     = 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        l0l1b_configuration(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

        l0l1b_config['settings']['van_cittert_steps'] = 0                 #switch off stray light
        l0l1b_config['filename']['level0']  ='Tango_Carbon_l0_exp7.2_1kernel_b1_'                   #level 1a input
        l0l1b_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_1kernel_not_corrected_b1_'
        l0l1b_config['filename']['ckd']     = 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        l0l1b_configuration(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

    # ======= L1 to L2 processor ===================================
    
    if(global_config['setup']['flag_l1bl2']):
        # choose baseline L1BL2 config
        shutil.copyfile(paths.project+paths.L1L2_module + 'l1bl2_config_baseline.yaml',
                        paths.project+paths.L1L2_module + 'l1bl2_config.yaml',)
        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config.yaml"))
        l1bl2_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_21kernel_corrected_b1_'
        l1bl2_config['filename']['level2'] = 'Tango_Carbon_l2_l1b_exp7.2_21kernel_corrected_b1_'
        l1bl2_config['pixel_mask'] = False
        l1bl2_config['isrf_settings']['type'] = 'Gaussian'                         #type of ISRF, currently only Gaussian or generalized_normal
        l1bl2_config['isrf_settings']['fwhm'] =  0.45                              #fwhm  [nm]

        level1b_to_level2_processor(paths, global_config, l1bl2_config)
 
        l1bl2_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_21kernel_not_corrected_b1_'
        l1bl2_config['filename']['level2']  = 'Tango_Carbon_l2_l1b_exp7.2_21kernel_not_corrected_b1_'
        l1bl2_config['pixel_mask'] = False
        level1b_to_level2_processor(paths, global_config, l1bl2_config)

        shutil.copyfile(paths.project+paths.L1L2_module + 'l1bl2_config_baseline.yaml',
                        paths.project+paths.L1L2_module + 'l1bl2_config.yaml',)
        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config.yaml"))
        l1bl2_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_1kernel_corrected_b1_'
        l1bl2_config['filename']['level2'] = 'Tango_Carbon_l2_l1b_exp7.2_1kernel_corrected_b1_'
        l1bl2_config['pixel_mask'] = False
        l1bl2_config['isrf_settings']['type'] = 'Gaussian'                         #type of ISRF, currently only Gaussian or generalized_normal
        l1bl2_config['isrf_settings']['fwhm'] =  0.45                              #fwhm  [nm]

        level1b_to_level2_processor(paths, global_config, l1bl2_config)
 
        l1bl2_config['filename']['level1b'] = 'Tango_Carbon_l1b_exp7.2_1kernel_not_corrected_b1_'
        l1bl2_config['filename']['level2']  = 'Tango_Carbon_l2_l1b_exp7.2_1kernel_not_corrected_b1_'
        l1bl2_config['pixel_mask'] = False
        level1b_to_level2_processor(paths, global_config, l1bl2_config)


    print('Experiment 7.2 sucessfully performed.')
    sys.exit()
