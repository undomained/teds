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
    temp_output_file = None
    # readin output file
    output_file = config['io'][output_key]
    output_data = dn.DataNetCDF(logger, output_file, mode='r')

    # cal_level determines until which step IM and L1B are run.
    # A check on cal_level can be used to determine whether or not reshaping is needed.

    # TODO: check if this works for both L1B and IM and possible in between steps
    if config['cal_level'] == 'rad':
        print("Not detector image. No need to reshape and create temporary output file")
        return temp_output_file
    elif config['cal_level'] == 'l1b':
        print("Not detector image. No need to reshape and create temporary output file")
        return temp_output_file
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
        temp_output_file = f"{output_file[:-3]}_temp.nc"
        print(f"temp_output: {temp_output_file}")
        output_data.write(temp_output_file)

    return temp_output_file

def get_file_name(config, step):
    """
        Determine output/input file path.
        When it is a nominal run it is just base_path
        When it concerns a scenarion and a step is required to be rerun 
        for this scenario it is base_path/scenarios/scenario_subdir
    """
    file_path = ''
    base_path = config['io']['base_dir']
    scenario_dir = ''
    if 'scenario' in config:
        scenario_dir = os.path.join('scenarios',config['scenario']['subdir'])
    output_path = os.path.join(base_path,scenario_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if ('scenario' in config) and (step in config['scenario']['steps']):
        file_path = output_path
    else:
        file_path = base_path
       
    return file_path 

def get_specific_config(logger, orig_config, kind):
    """
        Obtain module specific configuration from full configuration
    """
    if kind not in orig_config:
        error_msg = f"Unknown module {kind}. Unable to fetch corresponding config. Exiting!"
        logger.error(error_msg)
        sys.exit(error_msg)

    # Get the module specific settings
    specific_config = orig_config[kind]
    # Get the header
    specific_config['header'] = orig_config['header']
    # Get the module IO settings
    specific_config['io'] = {}

    if kind == 'GM':
        # Combine path and file name
        output_path = get_file_name(orig_config, 'gm')
        specific_config['io']['gm'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'SGM':
        # Combine path and file name
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm_rad'] = os.path.join(output_path, orig_config['io']['sgm_rad'])
        specific_config['io']['sgm_atm_raw'] = os.path.join(output_path, orig_config['io']['sgm_atm_raw'])
        specific_config['io']['sgm_atm'] = os.path.join(output_path, orig_config['io']['sgm_atm'])
        specific_config['io']['gm'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'L1L2':
        # Combine path and file name
        output_path = get_file_name(orig_config, 'l1l2')
        specific_config['io']['l2'] = os.path.join(output_path, orig_config['io']['l2'])

        output_path = get_file_name(orig_config, 'gm')
        specific_config['io']['gm'] = os.path.join(output_path, orig_config['io']['gm'])

        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm_atm'] = os.path.join(output_path, orig_config['io']['sgm_atm'])
        specific_config['io']['sgm_rad'] = os.path.join(output_path, orig_config['io']['sgm_rad'])

        output_path = get_file_name(orig_config, 'l1al1b')
        specific_config['io']['l1b'] = os.path.join(output_path, orig_config['io']['l1b'])

        # Also need acces to isrf which is a IM configuration parameter
        specific_config['isrf'] = orig_config['IM']['isrf']

    elif kind == 'L1AL1B':

        # Output to L1B
        output_path = get_file_name(orig_config, 'l1al1b')
        specific_config['io']['l1b'] = os.path.join(output_path, orig_config['io']['l1b'])

        # Input to L1B
        output_path = get_file_name(orig_config, 'im')
        l1a_file_name = orig_config['io']['l1a']

#        specific_config['io']['l1a'] = os.path.join(output_path, orig_config['io']['l1a'])

        # Note: if IM was run using the python code, the l1a python output needs to be readin
        do_python = orig_config['IM']['do_python']
        if do_python:
            l1a_file_name_python = f"{l1a_file_name[0:-4]}_python.nc"
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name_python)
        else:
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name)

        specific_config['io']['binning_table'] = orig_config['io']['binning_table']
        specific_config['io']['ckd'] = orig_config['io']['ckd']

    elif kind == 'IM':
        specific_config['io']['binning_table'] = orig_config['io']['binning_table']
        specific_config['io']['ckd'] = orig_config['io']['ckd_im']

        do_python = specific_config['do_python']
        # Output of IM
        output_path = get_file_name(orig_config, 'im')
        l1a_file_name = orig_config['io']['l1a']
        if do_python:
            l1a_file_name_python = f"{l1a_file_name[0:-3]}_python.nc"
            print(f"l1a python file name: {l1a_file_name_python}")
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name_python)
            specific_config['io']['im_algo_output'] = os.path.join(output_path, orig_config['io']['im_algo_output'])
        else:
#            specific_config['io']['l1a'] = os.path.join(output_path, orig_config['io']['l1a'])
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name)

        # Input to IM
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['l1b'] = os.path.join(output_path, orig_config['io']['l1b_im'])

    elif kind == 'PAM':
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm_atm'] = os.path.join(output_path, orig_config['io']['sgm_atm'])

        output_path = get_file_name(orig_config, 'l1l2')
        specific_config['io']['l2'] = os.path.join(output_path, orig_config['io']['l2'])

    else:
        error_msg = f"Unknown module kind {kind}. Unable to fetch corresponding IO config. Exiting!"
        logger.error(error_msg)
        sys.exit(error_msg)

    return specific_config

def add_module_specific_attributes(logger, config, attribute_dict, step):
    """
        Add the module specific settings to the attribute_dict
    """
    new_attribute_dict = attribute_dict.copy()

    for key, value in config.items():
        new_key = f"{step}_{key}"
        new_attribute_dict[new_key] = value

    return new_attribute_dict


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

        gm_config = get_specific_config(logger, configuration, 'GM')
        attribute_dict = add_module_specific_attributes(logger, gm_config, attribute_dict, 'gm')
        E2EModule = importlib.import_module("GM.gm")
        E2EModule.geometry_module(gm_config, logger=logger)
        # add attributes to the output file
        Utils.add_attributes_to_output(logger, gm_config['io']['gm'], attribute_dict)

    if step == 'sgm' or step == 'all':
        sgm_config = get_specific_config(logger, configuration, 'SGM')
        attribute_dict = add_module_specific_attributes(logger, sgm_config, attribute_dict, 'sgm')
        E2EModule = importlib.import_module("SGM.sgm_no2")
        E2EModule.scene_generation_module_nitro(logger,sgm_config)
        Utils.add_attributes_to_output(logger, sgm_config['io']['sgm_rad'], attribute_dict)
        Utils.add_attributes_to_output(logger, sgm_config['io']['sgm_atm'], attribute_dict)
        Utils.add_attributes_to_output(logger, sgm_config['io']['sgm_atm_raw'], attribute_dict)

    if step == 'im' or step == 'all':

        im_config = get_specific_config(logger, configuration, 'IM')
        attribute_dict = add_module_specific_attributes(logger, im_config, attribute_dict, 'im')
        # write config to temp yaml file with IM values filled in
        im_config_file = f"{cfg_path}/im_config_temp.yaml"
        with open(im_config_file,'w') as outfile:
            yaml.dump(im_config, outfile)

        if not im_config['do_python']:
            # Run C++ IM code
            # Need to call C++ using the IM specific config_file
            subprocess.run(["IM/tango_im/build/tango_im.x", im_config_file])

            # output dataset is 2D. In case of detector image (in case of some inbetween steps 
            # and of the final output) the second dimension is detector_pixels which is too large to view. 
            # Need to reshape to 3D to be able to make sense of this.
            temp_output_file = reshape_output(logger, 'l1a', im_config)

            # Add attributes to output file
            if temp_output_file is not None:
                Utils.add_attributes_to_output(logger, temp_output_file, attribute_dict)
            Utils.add_attributes_to_output(logger, im_config['io']['l1a'], attribute_dict)
        else:
            # run Python code
            E2EModule = importlib.import_module("IM.Python.instrument_model")
            E2EModule.instrument_model(im_config, logger, attribute_dict)
            # No need to reshape because Python output is already 3D

    if step == 'l1al1b' or step == 'all':

        l1b_config = get_specific_config(logger, configuration, 'L1AL1B')
        attribute_dict = add_module_specific_attributes(logger, l1b_config, attribute_dict, 'l1al1b')
        # write config to temp yaml file with L1B values filled in
        l1b_config_file = f"{cfg_path}/l1b_config_temp.yaml"
        with open(l1b_config_file,'w') as outfile:
            yaml.dump(l1b_config, outfile)

        # Need to call C++ with L1B specific config file
        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b.x", l1b_config_file])

        # output dataset is 2D. In case of detector image (in case of inbetween step)
        # the second dimension is detector_pixels which is too large to view. 
        # Need to reshape to 3D to be able to make sense of this.
        temp_output_file = reshape_output(logger, 'l1b', l1b_config)

        # Add attributes to output file
        if temp_output_file is not None:
            Utils.add_attributes_to_output(logger, temp_output_file, attribute_dict)
        Utils.add_attributes_to_output(logger, l1b_config['io']['l1b'], attribute_dict)

    if step == 'l1l2' or step == 'all':
        l2_config = get_specific_config(logger, configuration, 'L1L2')
        attribute_dict = add_module_specific_attributes(logger, l2_config, attribute_dict, 'l1l2')
        E2EModule = importlib.import_module("L1L2.l1bl2_no2")
        E2EModule.l1bl2_no2(logger, l2_config)
        Utils.add_attributes_to_output(logger, l2_config['l2'], attribute_dict)

    if step == 'pam' or step == 'all':
        pam_config = get_specific_config(logger, configuration, 'PAM')
        attribute_dict = add_module_specific_attributes(logger, pam_config, attribute_dict, 'pam')
        E2EModule = importlib.import_module("PAM.pam")
        E2EModule.pam_nitro(logger, pam_config)

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

