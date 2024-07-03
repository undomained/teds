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

def reshape_output(logger, output_key, config):
    """
        Note: output dataset is 2D. When it is on detector level
        these dimensions are: along track and pixels.
        This second dimension is too big (det_row x det_col)
        This makes viewing difficult.
        Read in the data and reshape to 3D in case of detector level
    """
    # readin output file
    output_file = config['io'][output_key]
    output_data = dn.DataNetCDF(logger, output_file, mode='r')

    # cal_level determines until which step IM and L1B are run.
    # A check on cal_level can be used to determine whether or not reshaping is needed.

    # TODO: check if this works for both L1B and IM and possible in between steps
    if config['cal_level'] == 'rad':
        print("Not detector image. No need to reshape and create temporary output file")
        return
    elif config['cal_level'] == 'l1b':
        print("Not detector image. No need to reshape and create temporary output file")
        return
    else:
        # Get data
        data = output_data.get('detector_image',  group='science_data', kind = 'variable')
        std_data = output_data.get('detector_stdev',  group='science_data', kind = 'variable')

        # Get dimensions from ckd?
        ckd_file = config['io']['ckd']
        ckd_data = dn.DataNetCDF(logger, ckd_file, mode='r')
        det_rows = ckd_data.get('detector_row', kind='dimension')
        det_cols = ckd_data.get('detector_column', kind='dimension')

        # At them to the output    
        output_data.add(name='detector_row', value=det_rows, kind='dimension')
        output_data.add(name='detector_column', value=det_cols, kind='dimension')
        output_data.add(name='along_track', value=data.shape[0], kind='dimension')
   
        # Data might be binned 
        # How do I get the right dimensions when binning has been applied?
        # For now stupidly devide det_rows by table_id
        bin_file = config['io']['binning_table']
        bin_data = dn.DataNetCDF(logger, bin_file, mode='r')
        bin_id = config['detector']['binning_table_id']
        table = f"Table_{bin_id}"
        print(f"Binning table: {table}")
        binned_pixels = bin_data.get('bins', group=table, kind='dimension')
        binned_rows = int(binned_pixels/det_cols)
        output_data.add(name='binned_row', value=binned_rows, kind='dimension')

        # is output level with or without binning? Check original shape of data
        n_pixels = data.shape[1]
        if int(n_pixels/det_cols) == det_rows:
            # detector image
            data_reshaped = np.reshape(data, (data.shape[0], det_rows, det_cols))
            output_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','detector_row','detector_column'), kind='variable')
            # standard deviation data is not always present. Check
            if std_data is not None:
                std_data_reshaped = np.reshape(std_data, (std_data.shape[0], det_rows, det_cols))
                output_data.add(name='detector_stdev_3d', value=std_data_reshaped, group='science_data', dimensions=('along_track','detector_row','detector_column'), kind='variable')
        elif int(n_pixels/det_cols) == binned_rows:
            # binned detector image
            data_reshaped = np.reshape(data, (data.shape[0], binned_rows, det_cols))
            output_data.add(name='detector_image_3d', value=data_reshaped, group='science_data', dimensions=('along_track','binned_row','detector_column'), kind='variable')
            # standard deviation data is not always present. Check
            if std_data is not None:
                std_data_reshaped = np.reshape(std_data, (data.shape[0], binned_rows, det_cols))
                output_data.add(name='detector_stdev_3d', value=std_data_reshaped, group='science_data', dimensions=('along_track','binned_row','detector_column'), kind='variable')

        print(f"SHAPE of Reshaped data image: {data_reshaped.shape}")
        # Remove the 2D images
        output_data.remove('detector_image', group='science_data', kind='variable')
        output_data.remove('detector_stdev', group='science_data', kind='variable')

        print("###################################################################")
        print("###################################################################")
        print(f"{output_data}")
        print("###################################################################")
        print("###################################################################")
        # write to temp output file
        temp_output = f"{output_file[:-3]}_temp.nc"
        print(f"temp_output: {temp_output}")
        output_data.write(temp_output)

    return

def create_im_and_l1b_specific_config_file(logger, configuration, specific_key):
    """
        A lot of settings for IM and L1B are the same.
        For some scenarios it is needed to be able to use different CKD file
        and different settings.
        For this purpose two keys(dicts) are defined: instrument_model and l1b
        Check which one is needed and move the settings one level up in the configuration.
        To avoid confusion remove the twe keys from config
        Write out to temporary config yaml file
    """

    # Copy the original configuration
    spec_configuration = configuration.copy()

    l1b_ckd = configuration['io']['ckd']
    l1b_file = configuration['io']['l1b']
    im_ckd = configuration['io']['ckd_im']
    im_input = configuration['io']['l1b_im']

    if specific_key == 'instrument_model':
        # Ensure we are using the IM CKD and l1b files (in IM it means the input file)
        spec_configuration['io']['ckd'] = im_ckd
        spec_configuration['io']['l1b'] = im_input
    elif specific_key == 'l1b':
        # Ensure we are using the L1B CKD and l1b files (in L1B it means the output file)
        spec_configuration['io']['ckd'] = l1b_ckd
        spec_configuration['io']['l1b'] = l1b_file
    else:
        error_msg = f"Specific key {specific_key} not understood. Exiting"
        logger.error(error_msg)
        sys.exit(error_msg)

    # Get the L1B or IM specific settings and move them one level up in configuration
    spec_config = configuration[specific_key]
    for key, value in spec_config.items():
        spec_configuration[key] = value

    # avoid confusion and delete IM and L1B specific settings from config file
    del spec_configuration['instrument_model']
    del spec_configuration['l1b']

    # write config to temp yaml file with IM or L1B values filled in
    spec_config_file = f"{cfg_path}/{specific_key}_config_temp.yaml"
    with open(spec_config_file,'w') as outfile:
        yaml.dump(spec_configuration, outfile)

    return spec_configuration, spec_config_file

def build(logger, config, step, cfg_path, attribute_dict):
    """
        Run E2E processor.
        - logger: Reference to the program logger
        - config: configuration file containing the settings for the different steps in the E2E processor
        - step: indicating which step in the E2E processor to run. 
        - cfg_path: configuration path
        - attribute_dict: Dictionary with attributes to be added to main of output netCDF files
    """

    # Make a copy of the full config file.
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

    if step == 'im' or step == 'all':

        # A lot of settings are shared between IM and L1B
        # For some scenarios it is needed to define them differently for IM and L1B
        # Therefore, two specific keys in config file. 
        # Move the settings corresponding to the instrument_model key one level up
        # in config file and write temp config out to file
        im_configuration, im_config_file = create_im_and_l1b_specific_config_file(logger, configuration, 'instrument_model')

        # Need to call C++ using the IM specific config_file
        subprocess.run(["IM/tango_im/build/tango_im.x", im_config_file])

        # output dataset is 2D. In case of detector image (in case of some inbetween steps 
        # and of the final output) the second dimension is detector_pixels which is too large to view. 
        # Need to reshape to 3D to be able to make sense of this.
        reshape_output(logger, 'l1a', im_configuration)

        Utils.add_attributes_to_output(logger, im_configuration['io']['l1a'], attribute_dict)


    if step == 'l1al1b' or step == 'all':

        # A lot of settings are shared between IM and L1B
        # For some scenarios it is needed to define them differently for IM and L1B
        # Therefore, two specific keys in config file. 
        # Move the settings corresponding to the l1b key one level up 
        # in config file and write temp config out to file
        l1b_configuration, l1b_config_file = create_im_and_l1b_specific_config_file(logger, configuration, 'l1b')

        # Need to call C++ with L1B specific config file
        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b.x", l1b_config_file])

        # output dataset is 2D. In case of detector image (in case of inbetween step)
        # the second dimension is detector_pixels which is too large to view. 
        # Need to reshape to 3D to be able to make sense of this.
        reshape_output(logger, 'l1b', l1b_configuration)

        # Add attributes to output file
        Utils.add_attributes_to_output(logger, l1b_configuration['io']['l1b'], attribute_dict)

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

