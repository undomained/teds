
# exp7.2: stray light analysis for a binary transition scene and with spectral dependence

# define  path to search for module

# define  path to search for module
import sys
path = "../"
if(str(path) not in sys.path):
    sys.path.append(path)

path = "../end_to_end/lib/"
if(str(path) not in sys.path):
    sys.path.append(path)
    
# import E2ES modules 
from end_to_end.lib import paths
from end_to_end.GM.gm import geometry_module
from end_to_end.SGM.sgm import scene_generation_module
from end_to_end.IM.create_im_configuration_file import im_configuration
from end_to_end.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
from end_to_end.L1L2.l1bl2 import level1b_to_level2_processor


import yaml
import sys
import subprocess

# ====================configuration part ======================================

class Emptyclass:
    pass

#run id
run_id = 'exp7.2'

# paths and file names
locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
locations.sgm['S2_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_S2_' + run_id + '.npy'
locations.sgm['afgl_input'] = paths.project + \
    paths.data_afgl + 'prof.AFGL.US.std'
locations.sgm['microHH_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_microHH_' + run_id + '.pkl'
locations.sgm['microHH_data_path'] = paths.project + \
    paths.data_microHH + 'Jaenschwalde_simulation1/'
locations.sgm['xsec_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_xsec_' + run_id + '.pkl'
locations.sgm['sun_reference'] = paths.project + paths.data_sol_spec + \
    'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['rad_output'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
locations.sgm['geo_output'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_atmosphere_' + run_id + '.nc'
locations.sgm['hapi_path'] =  paths.project + paths.data_harpi

locations.__setattr__('im', {})
locations.im['ckd_input'] = paths.project + paths.data_interface + \
    paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.im['rad_input'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
locations.im['binning_table'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'binning_table.nc'
locations.im['ckd_expo_time'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'optimal_exposure.txt'
locations.im['ckd_coadding'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'optimal_coadding.txt'
locations.im['IM_path'] = paths.project + paths.IM_module
locations.im['output'] = paths.project + paths.data_interface + \
    paths.interface_l1a + 'Tango_Carbon_l1a_' + run_id + '.nc'
locations.im['gm_input'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

locations.__setattr__('l1al1b', {})
locations.l1al1b['ckd_input']     = paths.project + paths.data_interface + \
    paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.l1al1b['binning_table'] = paths.project + paths.data_interface + \
    paths.interface_ckd + 'binning_table.nc'
locations.l1al1b['l1a_input']     = paths.project + paths.data_interface + \
    paths.interface_l1a + 'Tango_Carbon_l1a_' + run_id + '.nc'
locations.l1al1b['l1b_output']    = paths.project + paths.data_interface + \
    paths.interface_l1b + 'Tango_Carbon_l1b_' + run_id + '.nc'
locations.l1al1b['gm_input']      = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
locations.l1al1b['L1AL1B_path']   = paths.project + paths.L1AL1B_module

locations.__setattr__('l1bl2', {})
locations.l1bl2['l1b_input'] = paths.project + paths.data_interface + \
    paths.interface_l1b + 'Tango_Carbon_l1b_' + run_id + '.nc'
locations.l1bl2['pixel_mask'] = ''  #needs to be specified only if mask =True
locations.l1bl2['afgl_input'] = paths.project + \
    paths.data_afgl + 'prof.AFGL.US.std'
locations.l1bl2['sun_reference'] = paths.project + paths.data_sol_spec + \
    'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.l1bl2['l2_output'] = paths.project + paths.data_interface + \
    paths.interface_l2 + 'Tango_Carbon_l2_' + run_id + '.nc'
locations.l1bl2['xsec_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_xsec_l2_' + run_id + '.pkl'
locations.l1bl2['l2_diags'] = ''
locations.l1bl2['hapi_path'] =  paths.project + paths.data_harpi

#scene specification for profile single_pixel and swath
scene_spec = {}
scene_spec['numb_atm_scenes']= 2
scene_spec['scene_trans_index'] = [0, 50, 100]
scene_spec['sza'] = [70., 70.]
scene_spec['saa']= [0., 0.]
scene_spec['vza']= [0., 0.]
scene_spec['vaa']= [0., 0.]
scene_spec['albedo']= [0.15, 0.70]

# =============================================================================
#
#  Current version of the Tango E2E simulator can handle the following profiles
#   1.  individual_spectra: In this case only one single spectrum is processed.
#       The IM and L1 processor cannot be used in this case but are replaced by
#       the simplified moduel siml1. scene_spec required
#   2.  single_swath: This profile acounts for a single swath and a homogenous
#       reference scene. It is meant to test the IM and the L1 processor
#       scene_spec required
#   3.  S2_microHH: This simulates a full data granuale of 30 x 30 km2 using
#       flexibale geometry information and Sentinel 2 albedo data and
#       microHH CO2 plume simulations. scene_spec not required
# =============================================================================

profile = 'single_swath'  

settings= {}
settings['gm']    = True
settings['sgm']   = True
settings['im']    = True
settings['l1al1b']= True
settings['l1bl2'] = True
settings['save_yaml'] = True

if __name__ == "__main__":
    # end to end global config file
 
    # ======= geometry module ======================================
    # choose baseline GM config
    if(settings['gm']):

        config= yaml.safe_load(open(paths.project+paths.GM_module+'gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config, **scene_spec}
        gm_config['profile'] = profile
        if(settings['save_yaml']):
            filename = paths.project+paths.GM_module+'gm_config_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(gm_config,f,sort_keys=False) 

        geometry_module(gm_config)

    # ======= scene generator module ===============================
    # choose baseline SGM config
    if(settings['sgm']):
                
        sgm_config= yaml.safe_load(open(paths.project+paths.SGM_module + "sgm_config_baseline.yaml"))
        sgm_config = {**locations.sgm, **sgm_config, **scene_spec}
        sgm_config['profile'] = profile
        if(settings['save_yaml']):
            filename = paths.project+paths.SGM_module+'sgm_config_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(sgm_config,f,sort_keys=False) 

        scene_generation_module(sgm_config)

    # ======= The instrument model =================================        
    if(settings['im']):
        
        #with 21 stray light kernels
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id'] = 1
        im_config['settings']['sw_stray'] = 1
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_21kernel_b1_' + run_id + '.nc'
        im_configuration(locations.im, im_config)     

        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id'] = 1
        im_config['settings']['sw_stray'] = 1
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_1kernel_b1_' + run_id + '.nc'
        im_configuration(locations.im, im_config)     

        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #without straylight
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id'] = 1
        im_config['settings']['sw_stray'] = 0
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_wout_strayl_b1_' + run_id + '.nc'
        im_configuration(locations.im, im_config)     

        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)
        
    # ======= The L0L1 pocessor ====================================
    if(settings['l1al1b']):

        #stray light corrected , 21 kernels
        l1al1b_config = yaml.safe_load(open(paths.project+paths.L1AL1B_module + "l1al1b_config.yaml"))
        l1al1b_config['settings']['van_cittert_steps']   = 4 
        locations.l1al1b['ckd_input']     = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
        locations.l1al1b['l1a_input']     = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_21kernel_b1_' + run_id + '.nc'
        locations.l1al1b['l1b_output']    = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_21kernel_corrected_b1_' + run_id + '.nc'
        l1al1b_configuration(locations.l1al1b, l1al1b_config)     

        cmd_str = paths.project+paths.L1AL1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L1AL1B_module +'l1al1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L1AL1B_module)

        #stray light not corrected , 21 kernels
        l1al1b_config['settings']['van_cittert_steps'] = 0                 #switch off stray light
        locations.l1al1b['l1a_input']     = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_21kernel_b1_' + run_id + '.nc'
        locations.l1al1b['l1b_output']    = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_21kernel_not_corrected_b1_' + run_id + '.nc'
        l1al1b_configuration(locations.l1al1b, l1al1b_config)     

        cmd_str = paths.project+paths.L1AL1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L1AL1B_module +'l1al1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L1AL1B_module)

        #stray corrected , 1 kernels
        l1al1b_config = yaml.safe_load(open(paths.project+paths.L1AL1B_module + "l1al1b_config.yaml"))
        l1al1b_config['settings']['van_cittert_steps'] = 4
        locations.l1al1b['ckd_input']     = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.l1al1b['l1a_input']     = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_1kernel_b1_' + run_id + '.nc'
        locations.l1al1b['l1b_output']    = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_1kernel_corrected_b1_' + run_id + '.nc'
        l1al1b_configuration(locations.l1al1b, l1al1b_config)         

        cmd_str = paths.project+paths.L1AL1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L1AL1B_module +'l1al1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L1AL1B_module)

        #stray light not corrected , 1 kernels
        l1al1b_config = yaml.safe_load(open(paths.project+paths.L1AL1B_module + "l1al1b_config.yaml"))
        l1al1b_config['settings']['van_cittert_steps'] = 0                 #switch off stray light
        locations.l1al1b['ckd_input']     = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.l1al1b['l1a_input']     = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_1kernel_b1_' + run_id + '.nc'
        locations.l1al1b['l1b_output']    = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_1kernel_not_corrected_b1_' + run_id + '.nc'
        l1al1b_configuration(locations.l1al1b, l1al1b_config)         

        cmd_str = paths.project+paths.L1AL1B_module + 'tango_l1b/build/tango_l1b '+ paths.project+paths.L1AL1B_module +'l1al1b_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.L1AL1B_module)

    # ======= L1 to L2 processor ===================================
    
    if(settings['l1bl2']):
        # choose baseline L1BL2 config
        
        #stray light 21 kernel, corrected
        l1bl2_config= yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config_baseline.yaml"))
        l1bl2_config['pixel_mask']= False
        l1bl2_config['isrf_settings']['type']= 'Gaussian'  
        l1bl2_config['isrf_settings']['fwhm']=  0.45  
        l1bl2_config['profile'] = profile
        locations.l1bl2['l1b_input'] = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_21kernel_corrected_b1_' + run_id + '.nc'
        locations.l1bl2['l2_output'] = paths.project + paths.data_interface + \
            paths.interface_l2 + 'Tango_Carbon_l2_21kernel_corrected_b1_' + run_id + '.nc'
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}
        if(settings['save_yaml']):
            filename = paths.project+paths.L1L2_module+'l1l2_config_21kernel_corr_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(l1bl2_config,f,sort_keys=False) 

        level1b_to_level2_processor(l1bl2_config)

        #stray light 21 kernel, not corrected
        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config_baseline.yaml"))
        l1bl2_config['pixel_mask']= False
        l1bl2_config['isrf_settings']['type']= 'Gaussian'  
        l1bl2_config['isrf_settings']['fwhm']=  0.45  
        l1bl2_config['profile'] = profile
        locations.l1bl2['l1b_input'] = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_21kernel_not_corrected_b1_' + run_id + '.nc'
        locations.l1bl2['l2_output'] = paths.project + paths.data_interface + \
            paths.interface_l2 + 'Tango_Carbon_l2_21kernel_not_corrected_b1_' + run_id + '.nc'
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}
        if(settings['save_yaml']):
            filename = paths.project+paths.L1L2_module+'l1l2_config_21kernel_not_corr_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(l1bl2_config,f,sort_keys=False) 

        level1b_to_level2_processor(l1bl2_config)

        #stray light 1 kernel, corrected
        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config_baseline.yaml"))
        l1bl2_config['pixel_mask']= False
        l1bl2_config['isrf_settings']['type']= 'Gaussian'  
        l1bl2_config['isrf_settings']['fwhm']=  0.45  
        l1bl2_config['profile'] = profile
        locations.l1bl2['l1b_input'] = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_1kernel_corrected_b1_' + run_id + '.nc'
        locations.l1bl2['l2_output'] = paths.project + paths.data_interface + \
            paths.interface_l2 + 'Tango_Carbon_l2_1kernel_corrected_b1_' + run_id + '.nc'
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}
        if(settings['save_yaml']):
            filename = paths.project+paths.L1L2_module+'l1l2_config_1kernel_corr_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(l1bl2_config,f,sort_keys=False) 

        level1b_to_level2_processor(l1bl2_config)

        l1bl2_config = yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config_baseline.yaml"))
        l1bl2_config['pixel_mask']= False
        l1bl2_config['isrf_settings']['type']= 'Gaussian'  
        l1bl2_config['isrf_settings']['fwhm']=  0.45  
        l1bl2_config['profile'] = profile
        locations.l1bl2['l1b_input'] = paths.project + paths.data_interface + \
            paths.interface_l1b + 'Tango_Carbon_l1b_1kernel_not_corrected_b1_' + run_id + '.nc'
        locations.l1bl2['l2_output'] = paths.project + paths.data_interface + \
            paths.interface_l2 + 'Tango_Carbon_l2_1kernel_not_corrected_b1_' + run_id + '.nc'
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}
        if(settings['save_yaml']):
            filename = paths.project+paths.L1L2_module+'l1l2_config_1kernel_not_corr_'+run_id
            with open(f'{filename}.yaml', 'w',) as f :
                yaml.dump(l1bl2_config,f,sort_keys=False) 

        level1b_to_level2_processor(l1bl2_config)

    print('Experiment 7.2 sucessfully performed.')
    sys.exit()
