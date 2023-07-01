"""
This file executes a sensitivity study for a knowledge errors of the ISRF. 
We consxier two parameters of the generalized normal distribution function, 
    (1) parameter b describing the block shape of the ISRF for given FWHM
    (2) parameter a which is the FWHM
"""

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
from end_to_end.L1L2.l1bl2 import level1b_to_level2_processor
from end_to_end.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor

#import other modules
import yaml
import shutil

# ====================configuration part ======================================


class Emptyclass:
    pass

#run id
run_id = 'exp6'

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

locations.__setattr__('siml1b', {})
locations.siml1b['sgm_input']  = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
locations.siml1b['gm_input']   = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
locations.siml1b['l1b_output'] = paths.project + paths.data_interface + \
    paths.interface_l1b + 'Tango_Carbon_l1b_' + run_id + '.nc'

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
scene_spec['sza']    = [70., 60, 50, 40, 30, 20, 10, 0] 
scene_spec['saa']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene_spec['vza']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene_spec['vaa']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene_spec['albedo'] = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15] 

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

profile= 'individual_spectra'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']      = True
settings['sgm']     = True
settings['siml1b']  = True
settings['l1bl2']   = True
settings['sw_isrf_block'] = True
settings['sw_isrf_fwhm']  = True
# ====================main part ================================================
if __name__ == "__main__":

    # ======= geometry module ======================================

    if(settings['gm']):

        config= yaml.safe_load(open(paths.project+paths.GM_module+'gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config, **scene_spec}
        gm_config['profile'] = profile
        
        geometry_module(gm_config)

    # ======= scene generator module ===============================

    if(settings['sgm']):

        sgm_config= yaml.safe_load(open(paths.project+paths.SGM_module + "sgm_config_baseline.yaml"))
        sgm_config = {**locations.sgm, **sgm_config, **scene_spec}
        sgm_config['profile'] = profile

        scene_generation_module(sgm_config)

    # ======= The simplified IM and L1B model ======================
    if(settings['siml1b']):
        
        siml1b_config= yaml.safe_load(open(paths.project+paths.SIML1B_module + "siml1b_config_baseline.yaml"))
        siml1b_config = {**locations.siml1b, **siml1b_config}
        siml1b_config['profile'] = profile
        simplified_instrument_model_and_l1b_processor(siml1b_config)
    
    # ======= L1 to L2 processor ===================================
    if(settings['l1bl2']):
        # choose baseline L1BL2 config

        l1bl2_config= yaml.safe_load(open(paths.project+paths.L1L2_module + "l1bl2_config.yaml"))
        l1bl2_config['pixel_mask']= False
        l1bl2_config['isrf_settings']['type'] =  'generalized_normal' 
        l1bl2_config = {**locations.l1bl2, **l1bl2_config}

        if(settings['sw_isrf_block']):
            l1bl2_config['isrf_settings']['fwhm'] = 0.45
            for ibcoeff in range(0,11):
    
                if(ibcoeff == 0):
                    l1bl2_config['xsec_forced'] = True
                else:   
                    l1bl2_config['xsec_forced'] = False

                l1bl2_config['isrf_settings']['bcoeff'] = 0.4 + ibcoeff*0.01

                print('========================================')
                print('fwhm: ',l1bl2_config['isrf_settings']['fwhm'],' bcoeff: ', l1bl2_config['isrf_settings']['bcoeff'])
                print('========================================')
    
                str_bcoeff = "%.3f" % (l1bl2_config['isrf_settings']['bcoeff'])

                l1bl2_config['l2_output'] = paths.project + paths.data_interface + \
                    paths.interface_l2 + 'Tango_Carbon_l2__bcoeff'+str_bcoeff + run_id + '.nc'

                level1b_to_level2_processor(l1bl2_config)

        if(settings['sw_isrf_fwhm']):

            l1bl2_config['isrf_settings']['bcoeff'] = 0.45

            for iacoeff in range(0,11):
    
                l1bl2_config['isrf_settings']['fwhm'] = 0.43 + 0.004*iacoeff
                
                if(iacoeff == 0):
                    l1bl2_config['xsec_forced'] = True
                else:   
                    l1bl2_config['xsec_forced'] = False
                    
                print('========================================')
                print('fwhm: ',l1bl2_config['isrf_settings']['fwhm'],' bcoeff: ', l1bl2_config['isrf_settings']['bcoeff'])
                print('========================================')

                str_acoeff = "%.3f" % (l1bl2_config['isrf_settings']['fwhm'])
    
                l1bl2_config['l2_output'] = paths.project + paths.data_interface + \
                    paths.interface_l2 + 'Tango_Carbon_l2__acoeff'+str_acoeff + run_id + '.nc'
        
                level1b_to_level2_processor(l1bl2_config)
                