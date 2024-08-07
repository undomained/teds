#==============================================================================
#     using the local config to generate the l1al1b_config.cfg file
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
#==============================================================================
import sys
import yaml

def l1al1b_configuration(local_config):     

    lines = []
    #========================man group=========================================
    lines.append('[main]\n') 
    # Process name. In the present implementation this is always l1b.
    lines.append('process = l1b\n')
    #location of ckd file
    lines.append('ckd_file_in = ' + local_config['ckd_input']+'\n')
    #location of ckd binning table
    lines.append('binningtable_filename =' + local_config['binning_table']+'\n')
    # log file path
    lines.append('log_file_path = '+str(local_config['L1AL1B_log_path']) + '\n')
    # instrument calibration choice (either spexone or tango)
    lines.append('instrument_cal = '+str(local_config['L1AL1B_instrument_cal']) + '\n')
    #========================l1b group=========================================
    lines.append('[l1b]\n') 
    # Iteration limit for the stray light deconvolution procedure. Set to
    # 0 turn off the stray light correction.
    lines.append('stray_van_cittert_steps = '+str(local_config['settings_L1AL1B']['van_cittert_steps'])+'\n')
    # Whether to determine geolocation. Use 0 for now.
    lines.append('geolocation = '+str(local_config['settings_L1AL1B']['geolocation']) +'\n')
    #Location of the nL1A product (input detector images)
#    lines.append('l1a_files = '+ local_config['l1a_input']+'\n')
    lines.append('l1a_files = '+ local_config['l1a_file']+'\n')
    #Location of the level1b data
    lines.append('outputfile = ' + local_config['l1b_file']+'\n')
    #location of gm input
    lines.append('geometry_file = ' + local_config['gm_file']+'\n')
    # enable sub-module of level 0 to 1b processor (1 = switch on, 0 = switch off)
    # dark current
    lines.append('dark_apply = '+str(local_config['settings_L1AL1B']['sw_dark'])+'\n')
    # non-linearity
    lines.append('nonlin_apply = '+str(local_config['settings_L1AL1B']['sw_nonlin'])+'\n')
    # pixel response non-uniformity
    lines.append('prnu_apply = '+str(local_config['settings_L1AL1B']['sw_prnu'])+'\n')
    # radiometric calibration
    lines.append('rad_apply = '+str(local_config['settings_L1AL1B']['sw_rad'])+'\n')

    # write IM config file 
    new_config = open(local_config['L1AL1B_cfg_path']+'l1al1b_config.cfg','w')
    new_config.writelines(lines)
    new_config.close()

    return

if __name__ == '__main__' and __package__ is None:

    config = yaml.safe_load(open(sys.argv[1]))
    l1al1b_configuration(config)
