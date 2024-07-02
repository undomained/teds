import os, sys
import argparse
import logging
import importlib
import subprocess
import yaml
import numpy as np
import teds.lib.lib_utils as Utils
import lib.data_netcdf.data_netcdf as dn

def cmdline(arguments):
    """             
        Get the command line arguments
        - arguments: command line arguments
                
        return: 
        - cfgFile: the configuration file
    """         
                
    usage = """Run the TANGO E2E processor.
               The configuration file contains the settings for each step in the E2E processor."""

    cfgHelp = """The configuration file needed to run the E2E processor. Possible choices: Geomery
              (gm), Scene Generation (sgm), Instrument model (im), L1A to L1B (l1al1b), L1B to L2
              (l1l2) or all steps (all)."""
    parser = argparse.ArgumentParser(description= usage)
    parser.add_argument( "cfgFile", metavar="FILE", help=cfgHelp)
    parser.add_argument( "step", metavar='STEP', choices=['gm','sgm','im','l1al1b','l1l2','pam','all'],
                       help="The steps that the E2E processor has to run.")
                    
    args = parser.parse_args(arguments)
    cfgFile = args.cfgFile
    step = args.step
    
    return cfgFile, step

def build(logger, config, step, cfg_path, attribute_dict):
    """
        Run E2E processor.
        - logger: Reference to the program logger
        - config: configuration file containing the settings for the different steps in the E2E processor
        - step: indicating which step in the E2E processor to run. 
        - cfg_path: configuration path
        - attribute_dict: Dictionary with attributes to be added to main of output netCDF files
    """

    configuration = config.copy()

    if step == 'gm' or step == 'all':

        E2EModule = importlib.import_module("GM.gm")
        E2EModule.geometry_module(config, logger=logger)
        # add attributes to the output file
        Utils.add_attributes_to_output(logger, config['gm_file'], attribute_dict)

    if step == 'sgm' or step == 'all':
        E2EModule = importlib.import_module("SGM.sgm_no2")
        E2EModule.scene_generation_module_nitro(logger,config)
        Utils.add_attributes_to_output(logger, config['sgm_rad_file'], attribute_dict)
        Utils.add_attributes_to_output(logger, config['sgm_atm_file'], attribute_dict)

    

    l1b_ckd = config['io']['ckd']
    l1b_file = config['io']['l1b']
    im_ckd = config['io']['ckd_im']
    im_input = config['io']['l1b_im']

    if step == 'im' or step == 'all':
        # IM and L1B have setttings in common, but need to be flexible and be able to use other ckd file
        # or switch on/off other steps
        # IM files are in io: ckd_im and l1b_im
        # other settings are in key instrument_model
        # need the info one level up in config yaml file. 
        # write out temp yaml file

        im_configuration = configuration.copy()
        # Ensure we are using the IM CKD and l1b files (in IM it means the input file)
#        config['io']['l1b'] = im_input
#        config['io']['ckd'] = im_ckd
        im_configuration['io']['l1b'] = im_input
        im_configuration['io']['ckd'] = im_ckd
        # Get the IM specific steps and move them one level up
        im_config = config['instrument_model']
        for key, value in im_config.items():
            im_configuration[key] = value

        del im_configuration['instrument_model']
        del im_configuration['l1b']

        # write config to temp yaml file with IM values filled in
        im_config_file = f"{cfg_path}/im_config_temp.yaml"
        with open(im_config_file,'w') as outfile:
            yaml.dump(im_configuration, outfile)

        # Create cfg file to be used for IM executable

#        E2EModule = importlib.import_module("IM.create_im_configuration_file_nitro")
#        E2EModule.im_configuration(config)
#        subprocess.run(["IM/tango_ckd_model/build/ckdmodel", f"{cfg_path}/im_config.cfg"])
#        subprocess.run(["IM/tango_im/build/tango_im.x", cfgFile])

        # Need to call C++
        subprocess.run(["IM/tango_im/build/tango_im.x", im_config_file])
        # readin output file
        l1a_file = im_configuration['io']['l1a']
        l1a_data = dn.DataNetCDF(logger, l1a_file, mode='r')
#        if config['cal_level'] == 'swath' or config['cal_level'] == 'rad':
        if im_configuration['cal_level'] == 'rad':
            print("Not detector image")
        else:
            # Get data
            data = l1a_data.get('detector_image',  group='science_data', kind = 'variable')
            # Get dimensions from ckd?
            ckd_file = im_configuration['io']['ckd']
            ckd_data = dn.DataNetCDF(logger, ckd_file, mode='r')
            det_rows = ckd_data.get('detector_row', kind='dimension')
            det_cols = ckd_data.get('detector_column', kind='dimension')
    
            l1a_data.add(name='detector_row', value=det_rows, kind='dimension')
            l1a_data.add(name='detector_column', value=det_cols, kind='dimension')
            l1a_data.add(name='along_track', value=data.shape[0], kind='dimension')
    
            # reshape data
            # How do I get the right dimensions when binning has been applied?
            # For now stupidly devide det_rows by table_id
            bin_file = im_configuration['io']['binning_table']
            bin_data = dn.DataNetCDF(logger, bin_file, mode='r')
            bin_id = im_configuration['detector']['binning_table_id']
            table = f"Table_{bin_id}"
            print(f"Binning table: {table}")
            binned_pixels = bin_data.get('bins', group=table, kind='dimension')
            binned_rows = int(binned_pixels/det_cols)
            l1a_data.add(name='binned_row', value=binned_rows, kind='dimension')

#            data_reshaped = np.reshape(data, (data.shape[0], det_rows, det_cols))
            data_reshaped = np.reshape(data, (data.shape[0], binned_rows, det_cols))
            print(f"Reshaped data shape : {data_reshaped.shape}")
            print(f"Reshaped data image 1: {data_reshaped[1,:,:]}")
            if im_configuration['cal_level'] != 'l1a':
                l1a_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','detector_row','detector_column'), kind='variable')
            else:
                l1a_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','binned_row','detector_column'), kind='variable')
    
            print("###################################################################")
            print("###################################################################")
            print(f"{l1a_data}")
            print("###################################################################")
            print("###################################################################")
            # write to temp output file
            l1a_data.write('l1a_cpp_temp.nc')

        Utils.add_attributes_to_output(logger, im_configuration['io']['l1a'], attribute_dict)


    if step == 'l1al1b' or step == 'all':
#        # Create cfg file to be used for L1AL1B executable
#        E2EModule = importlib.import_module("L1AL1B.create_l1a1b_configuration_file_nitro")
#        E2EModule.l1al1b_configuration(config)
#        # Need to call C++
#        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b", f"{cfg_path}/l1al1b_config.cfg"])

        l1b_configuration = configuration.copy()

        # Ensure we are using the IM CKD and l1b files (in IM it means the input file)
        l1b_configuration['io']['ckd'] = l1b_ckd
        l1b_configuration['io']['l1b'] = l1b_file
        # Get the L1B specific steps and move them one level up
        l1b_config = configuration['l1b']
        for key, value in l1b_config.items():
            l1b_configuration[key] = value

        del l1b_configuration['instrument_model']
        del l1b_configuration['l1b']

        # write config to temp yaml file with IM values filled in
        l1b_config_file = f"{cfg_path}/l1b_config_temp.yaml"
        with open(l1b_config_file,'w') as outfile:
            yaml.dump(l1b_configuration, outfile)

        # Need to call C++
#        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b.x", cfgFile])
        result = subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b.x", l1b_config_file])
#        result = subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b.x", l1b_config_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#        print("STDOUT:")
#        print(result.stdout)
#        print("STDERR:")
#        print(result.stderr)

        # readin output file
        l1b_file = l1b_configuration['io']['l1b']
        l1b_data = dn.DataNetCDF(logger, l1b_file, mode='r')

        print("###################################################################")
        print(f"{l1b_data}")
        print("###################################################################")
        if l1b_configuration['cal_level'] == 'rad':
            print("Not detector image")
        elif l1b_configuration['cal_level'] == 'l1b':
            print("Not detector image")
        else:
            # Get data
            data = l1b_data.get('detector_image',  group='science_data', kind = 'variable')
            # Get dimensions from ckd?
            ckd_file = l1b_configuration['io']['ckd']
            ckd_data = dn.DataNetCDF(logger, ckd_file, mode='r')
            det_rows = ckd_data.get('detector_row', kind='dimension')
            det_cols = ckd_data.get('detector_column', kind='dimension')
    
            l1b_data.add(name='detector_row', value=det_rows, kind='dimension')
            l1b_data.add(name='detector_column', value=det_cols, kind='dimension')
            l1b_data.add(name='along_track', value=data.shape[0], kind='dimension')
    
            # reshape data
            # How do I get the right dimensions when binning has been applied?
            # For now stupidly devide det_rows by table_id
            bin_file = l1b_configuration['io']['binning_table']
            bin_data = dn.DataNetCDF(logger, bin_file, mode='r')
            bin_id = l1b_configuration['detector']['binning_table_id']
            table = f"Table_{bin_id}"
            print(f"Binning table: {table}")
            binned_pixels = bin_data.get('bins', group=table, kind='dimension')
            binned_rows = int(binned_pixels/det_cols)
            l1b_data.add(name='binned_row', value=binned_rows, kind='dimension')

#            data_reshaped = np.reshape(data, (data.shape[0], det_rows, det_cols))
            data_reshaped = np.reshape(data, (data.shape[0], binned_rows, det_cols))
            print(f"Reshaped data shape : {data_reshaped.shape}")
            print(f"Reshaped data image 1: {data_reshaped[1,:,:]}")
            if l1b_configuration['cal_level'] != 'l1b':
                l1b_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','binned_row','detector_column'), kind='variable')
            else:
                l1b_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','detector_row','detector_column'), kind='variable')
    
            print("###################################################################")
            print("###################################################################")
            print(f"{l1b_data}")
            print("###################################################################")
            print("###################################################################")
            # write to temp output file
            l1b_data.write('l1b_cpp_temp.nc')

#        Utils.add_attributes_to_output(logger, l1b_configuration['io']['l1b'], attribute_dict)

    if step == 'l1l2' or step == 'all':
        E2EModule = importlib.import_module("L1L2.l1bl2_no2")
        E2EModule.l1bl2_no2(logger, config)
        Utils.add_attributes_to_output(logger, config['l2_file'], attribute_dict)

    if step == 'pam' or step == 'all':
        E2EModule = importlib.import_module("PAM.pam")
        E2EModule.pam_nitro(logger, config)

if __name__ == "__main__":

    
    # Get logger for run script
    build_logger = Utils.get_logger()

    # Get input arguments
    cfgFile, step =  cmdline(sys.argv[1:])
    cfg_path, filename = os.path.split(cfgFile)

    # Get configuration info
    config = Utils.getConfig(build_logger, cfgFile)
    config['header']['path'] = cfg_path

    # Get information (like git hash and config file name and version (if available) 
    # that will be added to the output files as attributes
    main_attribute_dict = Utils.get_main_attributes(config)

    build(build_logger, config, step, cfg_path, main_attribute_dict)

