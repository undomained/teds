#==============================================================================
# using the config dictionary to generate a yaml file with comments
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys

def create_sgm_config_file(filename, config):     

    lines = []
    #==========================main============================================
    
    lines.append('# gm and sgm profile \n') 
    lines.append('profile:                 ' + config['profile'] + '  # must be the same for gm and sgm\n')
    lines.append('# file names and paths \n') 
    lines.append('# sgm output file \n') 
    lines.append('gm_input:                ' + config['gm_input'] + '\n')
    lines.append('afgl_input:              ' + config['afgl_input'] + '\n')
    lines.append('rad_putput:              ' + config['rad_output'] + '\n')
    lines.append('geo_output:              ' + config['geo_output'] + '\n')
    lines.append('sun_reference:           ' + config['sun_reference'] + '\n')
    lines.append('S2_dump:                 ' + config['S2_dump'] + '\n')
    lines.append('xsec_dump:               ' + config['xsec_dump'] + '\n')
    lines.append('microHH_dump:            ' + config['microHH_dump'] + '\n')
    lines.append('microHH_data_path:       ' + config['microHH_data_path'] + '\n')
    lines.append('hapi_path:               ' + config['hapi_path'] + '\n')
    lines.append('\n')


    if(config['profile']=='S2_microHH'):
        
        lines.append('# define 2D kernel for spatial smoothing og the data  \n') 
        lines.append('kernel_parameter: ' + '\n')
        lines.append('  type:              ' + str(config['kernel_parameter']['type']) + '  # type pf kernel, current version has only one option  \n')
        lines.append('  fwhm_x:            ' + str(config['kernel_parameter']['fwhm_x']) + '  # FWHM in ACT direction [m]\n')
        lines.append('  fwhm_y:            ' + str(config['kernel_parameter']['fwhm_y']) + '  # FWHM in ALT direction [m]\n')
        lines.append('  size_factor:       ' + str(config['kernel_parameter']['size_factor']) + '  # the size of the convolution kernel in units of FWHM \n')
        lines.append('\n')

        lines.append('# vertical grid parameters of the model atmosphere  \n') 
        lines.append('atmosphere: \n')
        lines.append('  nlay:              ' + str(config['atmosphere']['nlay']) + '  # number of atmospheric layers \n')
        lines.append('  dzlay:             ' + str(config['atmosphere']['dzlay']) + '  # geometrical thickness of the layers [m] \n')
        lines.append('\n')
                
        lines.append('# AFGL model atmosphere settings  \n') 
        lines.append('only_afgl:           ' + str(config['only_afgl']) + '  # switch to ignore microHH and use the AFGL atmosphere over the entire domain \n')
        lines.append('\n')
                
        lines.append('# microHH settings  \n') 
        lines.append('microHH: ' + '\n')
        lines.append('  time_stamp:        ' + str(config['microHH']['time_stamp']) + '  # time stamp of microHH simulation \n')
        lines.append('  lat_lon_src:       ' + str(config['microHH']['lat_lon_src']) + '  # type of microHH simulation for given target area \n')                
        lines.append('\n')

        lines.append('# spectral settings  \n') 
        lines.append('spec_settings: ' + '  # settings for line-by-line spectra \n')
        lines.append('  wave_start:        ' + str(config['spec_settings']['wave_start']) + '  # start wavelength [nm] \n')
        lines.append('  wave_end:          ' + str(config['spec_settings']['wave_end']) + '  # end wavelength [nm] \n')
        lines.append('  dwave:             ' + str(config['spec_settings']['dwave']) + '  # sampling distance [nm] \n')
        lines.append('\n')

        lines.append('# switches to force xsec, S2, and microHH simulation  \n') 
        lines.append('xsec_forced:         ' + str(config['xsec_forced']) + '  # flag to force new cross section calculation\n')
        lines.append('s2_forced:           ' + str(config['s2_forced']) + '  # flag to force Sentinel 2 calculation\n')
        lines.append('microHH_forced:      ' + str(config['microHH_forced']) + '  # flag to force microHH data processing\n')
        lines.append('\n')
                            
    if(config['profile']=='single_swath'):
        lines.append('# albedo variation within the scene, needs to be consistent with gm config \n') 
        lines.append('scene_spec: ' + '\n')
        lines.append('  numb_atm:          ' + str(config['scene_spec']['numb_atm']) + '  # number of different model atmospheres \n')
        lines.append('  scene_trans_index: ' + str(config['scene_spec']['scene_trans_index']) + '  # act index where the model atmosphere changes, should start and end with smallest and largest act index\n')
        lines.append('  albedo:            ' + str(config['scene_spec']['albedo']) + '  # Lambertian surface albedo [1] \n')

    # Location of SGM output (input spectra

    # write IM config file 
    sgm_config = open(filename,'w')
    sgm_config.writelines(lines)
    sgm_config.close()
    return
