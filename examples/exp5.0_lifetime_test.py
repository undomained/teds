
# exp4: single image performance with full IM and L0L1B
#       can be used for ckd degradation study

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

    dose = ['0.0', '1.8', '6.4', '8.6', '13.1', '15.5'] #different radiation doses in kRad
    
    global_config = yaml.safe_load(open("./exp5.0_lifetime_test.yaml"))

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
        ckd_file = 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['filename']['ckd']        = ckd_file
        im_config['noise']['switch']        = 0
        im_config['settings']['bin_id']     = 1
        im_config['settings']['sw_stray']   = 0
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    # ======= The L0L1 pocessor ====================================
    if(global_config['setup']['flag_l0l1b']):
        ckd_file = 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
        l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
        l0l1b_config['settings']['van_cittert_steps'] = 0  #
        l0l1b_config['filename']['ckd'] = ckd_file
        l0l1b_config['filename']['level1b'] = l0l1b_config['filename']['level1b']+'ref_'
        l0l1b_configuration(paths, global_config, l0l1b_config)     
        cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)

    for ds in dose:        #loop over differnet ckd files

        # ======= The instrument model =================================   
        if(global_config['setup']['flag_im']):
            ckd_file = 'OWL640S_low-gain_radiation_ckd_dose'+ds+'_21kernel'   
            im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
            im_config['filename']['ckd']        = ckd_file
            im_config['noise']['switch']        = 0
            im_config['settings']['bin_id']     = 1
            im_config['settings']['sw_stray']   = 0
            im_configuration(paths,global_config, im_config)     
            cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
            subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        
        # ======= The L0L1 pocessor ====================================
        if(global_config['setup']['flag_l0l1b']):
            ckd_file = 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
            l0l1b_config = yaml.safe_load(open(paths.project+paths.L0L1B_module + "l0l1b_config.yaml"))
            l0l1b_config['settings']['van_cittert_steps'] = 0  #
            l0l1b_config['filename']['ckd'] = ckd_file
            l0l1b_config['filename']['level0'] = im_config['filename']['level1a']
            l0l1b_configuration(paths, global_config, l0l1b_config)     
            cmd_str = paths.project+paths.L0L1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L0L1B_module +'l0l1b_config.cfg'
            subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L0L1B_module)
    
            file1b = l0l1b_config['filename']['level1b']+global_config['profile']+'_'+global_config['run_id']+'.nc'
            file1b_ref = l0l1b_config['filename']['level1b']+'ref_'+global_config['profile']+'_'+global_config['run_id']+'.nc'
            print(file1b_ref)
            l0l1b_config['filename']['pixel_mask'] = l0l1b_config['filename']['pixel_mask'] + 'T3_'+'dose_'+ds
            determine_pixel_mask(paths,l0l1b_config['filename']['pixel_mask'], file1b, file1b_ref)


        # ======= L1 to L2 processor ===================================
        if(global_config['setup']['flag_l1bl2']):
            # choose baseline L1BL2 config
            shutil.copyfile(paths.project+paths.L1L2_module + 'l1bl2_config_baseline.yaml',
                            paths.project+paths.L1L2_module + 'l1bl2_config.yaml',)
            l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config.yaml"))
            l1bl2_config['filename']['level2'] = 'Tango_Carbon_l2_dose'+ds
            l1bl2_config['pixel_mask'] = True
            l1bl2_config['filename']['pixel_mask'] = l1bl2_config['filename']['pixel_mask'] +'T3_'+ 'dose_'+ds
            l1bl2_config['isrf_settings']['type'] =  'Gaussian' 
            l1bl2_config['isrf_settings']['fwhm'] = 0.45

            level1b_to_level2_processor(paths, global_config, l1bl2_config)
 
    print('Experiment 5 sucessfully performed.')
    sys.exit()
