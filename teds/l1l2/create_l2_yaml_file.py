#==============================================================================
# using the config dictionary to generate a yaml file with comments
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys

def create_l1bl2_config_file(filename, config):     

    lines = []
    #==========================main============================================
    
    lines.append('# l1bl2 in/output file \n') 
    lines.append('l1b_input:          ' + config['l1b_input'] + '\n')
    lines.append('afgl_input:         ' + config['afgl_input'] + '\n')
    if(config['retrieval_init']['sw_pixel_mask']):
        lines.append('pixel_mask:         ' + config['pixel_mask'] + '\n')
    lines.append('sun_reference:      ' + config['sun_reference'] + '\n')
    lines.append('l2_output:          ' + config['l2_output'] + '\n')
    lines.append('xsec_dump:          ' + config['xsec_dump'] + '\n')
#    lines.append('l2_diags:           ' + config['l2_diags'] + '\n')
    lines.append('hapi_path:          ' + config['hapi_path'] + '\n')
    lines.append('\n')
    lines.append('retrieval_init:         #initialization parameter for inversion \n')
    lines.append('  max_iter:         ' + str(config['retrieval_init']['max_iter']) + '  # maximum number of iterations \n')
    lines.append('  chi2_lim:         ' + str(config['retrieval_init']['chi2_lim']) + '  # chi2 threshold value \n')
    lines.append('  sw_pixel_mask:    ' + str(config['retrieval_init']['sw_pixel_mask']) + '  # pixel mask flag\n')
    lines.append('  sw_ALT_select:    ' + str(config['retrieval_init']['sw_ALT_select']) + '  # flag for ALT image selection\n')
    if(config['retrieval_init']['sw_ALT_select']):
        lines.append('  first_ALT_index:  ' + str(config['retrieval_init']['first_ALT_index']) + '  # first ALT image index\n')
        lines.append('  last_ALT_index:   ' + str(config['retrieval_init']['last_ALT_index']) + '  # last ALT image index\n')
    lines.append('\n')
    lines.append('atmosphere:             #setting for model atmosphere \n') 
    lines.append('  nlay:             ' + str(config['atmosphere']['nlay'])  + '  # number of atmospheric layers\n')
    lines.append('  dzlay:            ' + str(config['atmosphere']['dzlay'])  + '  # geometrical thickness of the layers [m]\n')
    lines.append('  psurf:            ' + str(config['atmosphere']['psurf'])  + '  # surface pressure [Pa]\n')
    lines.append('\n')
    lines.append('spec_settings:              #spectra fit window \n')
    lines.append('  wavestart:        ' + str(config['spec_settings']['wavestart'])   +  ' # initial wavelength (indicator)\n')
    lines.append('  waveend:          ' + str(config['spec_settings']['waveend'])     +  ' # final wavlength of measurement\n')
    lines.append('  wave_extend:      ' + str(config['spec_settings']['wave_extend']) +  ' # shortwave and longwave extension of the measurement grid [nm]\n')
    lines.append('  dwave:            ' + str(config['spec_settings']['dwave'])       +  ' # spectral sampling of line-by-line spectra [nm]\n')
    lines.append('\n')
    lines.append('isrf_settings:          #settings for the instrument spectral response function \n')
    lines.append('  type:             ' + config['isrf_settings']['type'] + '  # type of ISRF, currently only Gaussian or generalized_normal\n')
    lines.append('  fwhm:             ' + str(config['isrf_settings']['fwhm']) + '  # fwhm / acoeff  [nm]\n')
    if(config['isrf_settings']['type']=='generalized_normal'):
        lines.append('  bcoeff:          ' + str(config['isrf_settings']['bcoeff']) + '  # only for generalized_normal, bcoeff=1/2 means Gaussian\n')
    lines.append('\n')
    lines.append('surface: \n')     
    lines.append('  wave_alb:         ' + str(config['surface']['wave_alb']) + '  # wavelengths at which the surface albedo is specified in the output\n')    
    lines.append('\n')
    lines.append('xsec_forced:        ' + str(config['xsec_forced']) + ' # flag to froce new cross section calculation\n')

    # write l1bl2 config file 
    l1bl2_config = open(filename,'w')
    l1bl2_config.writelines(lines)
    l1bl2_config.close()
    return
