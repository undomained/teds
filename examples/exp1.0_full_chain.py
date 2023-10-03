"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
# path = "../"
# if(str(path) not in sys.path):
#     sys.path.append(path)

# path = "../end_to_end/lib/"
# if(str(path) not in sys.path):
#     sys.path.append(path)
    
# import E2ES modules 
from end_to_end.lib import paths
from end_to_end.GM.gm import geometry_module
from end_to_end.GM.create_gm_yaml_file import create_gm_config_file
from end_to_end.SGM.sgm import scene_generation_module
from end_to_end.SGM.create_sgm_yaml_file import create_sgm_config_file
from end_to_end.IM.create_im_configuration_file import im_configuration
from end_to_end.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
from end_to_end.L1L2.l1bl2 import level1b_to_level2_processor
from end_to_end.L1L2.create_l2_yaml_file import create_l1bl2_config_file
from end_to_end.L1L2.SyntheticLevel2 import simplified_level2

#import other modules
import yaml
import subprocess

# ====================configuration part ======================================

class Emptyclass:
    pass

#run id
run_id = 'exp1.0_sgm_with_noise'

# paths and file names
locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
locations.sgm['S2_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_S2_exp1.0.npy'
#    paths.data_tmp + 'Tango_Carbon_S2_' + run_id + '.npy'
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
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_exp1.0_bin1.nc'
#    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
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
    paths.interface_gm + 'Tango_Carbon_gm_exp1.0_bin1.nc'
#    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

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
    paths.interface_gm + 'Tango_Carbon_gm_exp1.0_bin1.nc'
#    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
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
locations.l1bl2['sgm_input'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_atmosphere_' + run_id + '.nc'

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

profile= 'orbit'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']        = False
settings['sgm']       = False
settings['im']        = False
settings['l1al1b']    = False
settings['l1bl2']     = False
settings['save_yaml'] = False
settings['sl2']       = True
# ====================main part ================================================
if __name__ == "__main__":

    # ======= geometry module ======================================

    if(settings['gm']):

        config= yaml.safe_load(open(paths.project+paths.GM_module+'gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config}
        gm_config['profile'] = profile
        if(settings['save_yaml']):
            gm_yaml = paths.project+paths.GM_module+'gm_config_'+run_id + '.yaml'
            create_gm_config_file(gm_yaml, gm_config)
        geometry_module(gm_config)

    # ======= scene generator module ===============================

    if(settings['sgm']):

        sgm_config= yaml.safe_load(open(paths.project+paths.SGM_module + "sgm_config_baseline.yaml"))
        sgm_config = {**locations.sgm, **sgm_config}
        sgm_config['profile'] = profile
        if(settings['save_yaml']):
            sgm_yaml = paths.project+paths.SGM_module+'sgm_config_'+run_id + '.yaml'
            create_sgm_config_file(sgm_yaml, sgm_config)
        scene_generation_module(sgm_config)

    # ======= The instrument model =================================

    if(settings['im']):

        # with constant stray light kernel
        im_config= yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch']= 1
        im_config['settings']['bin_id']= 1
        im_config['settings']['sw_stray']= 0
        im_configuration(locations.im, im_config)
        cmd_str= paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    # ======= The L0L1 pocessor ====================================
    if(settings['l1al1b']):

        l1al1b_config= yaml.safe_load(open(paths.project+paths.L1AL1B_module + "l1al1b_config.yaml"))
        l1al1b_config['settings']['van_cittert_steps']= 0
        l1al1b_configuration(locations.l1al1b, l1al1b_config)
        cmd_str= paths.project+paths.L1AL1B_module + 'tango_l1b/build/tango_l1b ' + paths.project+paths.L1AL1B_module + 'l1al1b_config.cfg'
        subprocess.run(cmd_str, shell=True,
                       cwd=paths.project+paths.L1AL1B_module)

    # ======= L1 to L2 processor ===================================
    if(settings['l1bl2']):
        # choose baseline L1BL2 config

        l1bl2_config= yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config_baseline.yaml"))
        l1bl2_config['retrieval_init']['sw_pixel_mask']= False
        l1bl2_config['isrf_settings']['type']= 'Gaussian'  # type of ISRF, currently only Gaussian or generalized_normal
        l1bl2_config['isrf_settings']['fwhm']=  0.45  # fwhm  [nm]
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}
        if(settings['save_yaml']):
            l1bl2_yaml = paths.project+paths.L1L2_module+'l1bl2_config_'+run_id + '.yaml'
            create_l1bl2_config_file(l1bl2_yaml, l1bl2_config)
        level1b_to_level2_processor(l1bl2_config)

    # ======= L1 to L2 processor ===================================
    if(settings['sl2']):

        sl2_config = {}
        sl2_config['instrument'] = 'TANGO_Carbon'
        sl2_config['sgm_file'] = paths.project + paths.data_interface + paths.interface_sgm + 'Tango_Carbon_sgm_atmosphere_' + run_id + '.nc'
        sl2_config['l2_file'] = paths.project + paths.data_interface + paths.interface_l2 + 'Tango_Carbon_l2_' + run_id + '.nc'
        sl2_config['l2_file'] = paths.project + paths.data_interface + paths.interface_l2 + 'Tango_Carbon_l2_test.nc'
        sl2_config['precision relative'] = 0.005
        sl2_config['precision constant'] = 0.00
        sl2_config['seed'] = 10
        
        simplified_level2(sl2_config['sgm_file'], sl2_config['l2_file'], sl2_config['instrument'], sl2_config['precision relative'], sl2_config['precision constant'], sl2_config['seed'])
