#==============================================================================
# using the config dictionary to generate a yaml file with comments
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys

def create_gm_config_file(filename, config):     

    lines = []
    #==========================main============================================
    
    lines.append('# gm and sgm profile \n') 
    lines.append('profile:   ' + config['profile'] + '  # must be the same for gm and sgm\n')
    lines.append('# file names and paths \n') 
    lines.append('# gm output file \n') 
    lines.append('output:    ' + config['output'] + '\n')
    lines.append('\n')
    if(config['profile']=='S2_microHH'):
        lines.append('# define the field of regard in terms of number of pixels and ACT geometry  \n') 
        lines.append('field_of_regard: ' + '\n')
        lines.append('  nact:              ' +str(config['field_of_regard']['nact']) + '  # number of pixels in across track direction \n')
        lines.append('  nalt:              ' + str(config['field_of_regard']['nalt']) + '  # number of pixels in along track direction \n')
        lines.append('  alpha_act_min:     ' + str(config['field_of_regard']['alpha_act_min']) + '  # minimum ACT angle at instrument level \n')
        lines.append('  alpha_act_max:     ' + str(config['field_of_regard']['alpha_act_max']) + '  # maximum ACT angle at instrument level \n')
        lines.append('# define central point of data granule (lat,lon) \n') 
        lines.append('\n')
        lines.append('geometry: ' + '\n')
        lines.append('  lat_initial:       ' + str(config['geometry']['lat_initial']) + '  # latitude coordinate of central point \n')
        lines.append('  lon_initial:       ' + str(config['geometry']['lon_initial']) + '  # longitude coordinate of central point \n')
        lines.append('# time specifications \n') 
        lines.append('\n')
        lines.append('time: ' + '\n')
        lines.append('  time_incremental:  ' + str(config['time']['time_incremental']) + '  # time incremental for consecutive readouts for FMC = 1 \n')
        lines.append('  year:              ' + str(config['time']['year'])     + '  # date information to determine solar geometry \n')
        lines.append('  month:             ' + str(config['time']['month'])    +'\n')
        lines.append('  day:               ' + str(config['time']['day'])      +'\n')
        lines.append('  hour:              ' + str(config['time']['hour'])     +'\n')
        lines.append('  minute:            ' + str(config['time']['minute'])   +'\n') 
        lines.append('  timezone:          ' + str(config['time']['timezone']) +'\n') 
        lines.append('\n')
        lines.append('# satellite and instrument attitude data \n') 
        lines.append('satellite: ' + '\n')
        lines.append('  alpha_roll: ' + str(config['satellite']['alpha_roll']) + '  # roll angle [degree] \n')
        lines.append('  sat_height:        ' + str(config['satellite']['sat_height']) + '  # satellite height [m] \n')
        lines.append('\n')

    if(config['profile']=='single_swath'):
        lines.append('# setup up geometry and albedo of the scene \n') 
        lines.append('scene_spec: ' + '\n')
        lines.append('  numb_atm:          ' + str(config['scene_spec']['numb_atm']) + '  # number of different model atmospheres \n')
        lines.append('  scene_trans_index: ' + str(config['scene_spec']['scene_trans_index']) + '  # act index where the model atmosphere changes, should start and end with smallest and largest act index\n')
        lines.append('  sza:               ' + str(config['scene_spec']['sza']) + '  # solar zenith angle [degeee] \n')
        lines.append('  saa:               ' + str(config['scene_spec']['sza']) + '  # solar azimuth angle [degeee] \n')
        lines.append('  vza:               ' + str(config['scene_spec']['vza']) + '  # viewing zenith angle [degeee] \n')
        lines.append('  vaa:               ' + str(config['scene_spec']['vza']) + '  # viewing azimuth angle [degeee] \n')
        lines.append('  sza:               ' + str(config['scene_spec']['sza']) + '  # solar zenith angle [degeee] \n')
        lines.append('  saa:               ' + str(config['scene_spec']['sza']) + '  # solar azimuth angle [degeee] \n')
        lines.append('\n')

    if(config['profile']=='individual_spectra'):
        lines.append('# setup up geometry and albedo of the scene \n') 
        lines.append('  sza:               ' + str(config['scene_spec']['sza']) + '  # solar zenith angle [degeee] \n')
        lines.append('  saa:               ' + str(config['scene_spec']['sza']) + '  # solar azimuth angle [degeee] \n')
        lines.append('  vza:               ' + str(config['scene_spec']['vza']) + '  # viewing zenith angle [degeee] \n')
        lines.append('  vaa:               ' + str(config['scene_spec']['vza']) + '  # viewing azimuth angle [degeee] \n')
        lines.append('  sza:               ' + str(config['scene_spec']['sza']) + '  # solar zenith angle [degeee] \n')
        lines.append('  saa:               ' + str(config['scene_spec']['sza']) + '  # solar azimuth angle [degeee] \n')

    # Location of SGM output (input spectra

    # write IM config file 
    gm_config = open(filename,'w')
    gm_config.writelines(lines)
    gm_config.close()
    return
