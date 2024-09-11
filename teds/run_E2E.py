# pylint: disable=unused-import, fixme, logging-fstring-interpolation
"""Providing code to run the Tango E2E processor"""

import argparse
from datetime import timedelta
import importlib
import os
import subprocess
import sys
import time
import numpy as np
import yaml

from teds import log
import teds.lib.lib_utils as Utils
import teds.lib.data_netcdf.data_netcdf as dn
import logging as _logging

def cmdline(arguments):
    """
        Get the command line arguments

        : param arguments: command line arguments
        : type arguments: List

        :return: cfg_file: the configuration file path
        :rtype:  cfg_file: String
        :return: step: the step of the E2E which needs to be run
        :rtype:  step: String
    """
    usage = """Run the TANGO E2E processor.
               The configuration file contains the settings for each step in the E2E processor."""

    cfg_help = """The configuration file needed to run the E2E processor. Possible choices: Geomery
              (gm), Scene Generation (sgm), Instrument model (im), L1A to L1B (l1al1b), L1B to L2
              (l1l2) or all steps (all)."""
    parser = argparse.ArgumentParser(description= usage)
    parser.add_argument( "cfg_file", metavar="FILE", help=cfg_help)
    parser.add_argument( "step", metavar='STEP',
                         choices=['gm','sgm','im','l1al1b','l1l2','pam','all'],
                         help="The steps that the E2E processor has to run.")

    args = parser.parse_args(arguments)
    cfg_file = args.cfg_file
    step = args.step

    return cfg_file, step

def is_reshape_needed(output_key, config):
    """
        Check if reshaping of the data from 2D to 3D is needed.

        :param output_key: key to use to find output file in configuration
        :type output_key: String
        :param config: configuration
        :type config: Dictionary

        :return reshape_needed: to indicated if reshaping is needed
        :rtype reshape_needed: Boolean
    """
    reshape_needed = True
    if config['cal_level'] == 'rad':
        log.info("Radiometric level. Not detector image. No need to reshape and create temporary output file")
        reshape_needed = False
    elif output_key=='l1b' and config['cal_level'] == 'swath':
        log.info("L1B and Swath. Still detector image. No need to reshape and create temporary output file")
        reshape_needed = False
    elif config['cal_level'] == 'l1b':
        log.info("L1B level. Not detector image. No need to reshape and create temporary output file")
        reshape_needed = False

    return reshape_needed

def detector_image_dimensions(output_data, data, config ):
    """
        :param output_data: netcdf data
        :type output_data: dataNetcdf
        :param config: configuration
        :type config: Dictionary

        :return : List of detector row, detector columns and binned_rows
        :rtype  : List
    """

    # Get dimensions from ckd?
    ckd_file = config['io']['ckd']
    ckd_data = dn.DataNetCDF(ckd_file, mode='r')
    det_rows = ckd_data.get('detector_row', kind='dimension')
    det_cols = ckd_data.get('detector_column', kind='dimension')

    # At them to the output
    output_data.add(name='detector_row', value=det_rows, kind='dimension')
    output_data.add(name='detector_column', value=det_cols, kind='dimension')
# along_track dimension is alread part of the dataset. No need to add it.
#    output_data.add(name='along_track', value=data.shape[0], kind='dimension')

    # Data might be binned
    # How do I get the right dimensions when binning has been applied?
    # For now stupidly devide det_rows by table_id
    bin_file = config['io']['binning_table']
    bin_data = dn.DataNetCDF(bin_file, mode='r')
    if 'detector' in config:
        bin_id = config['detector']['binning_table_id']
    else:
        binning_ids = output_data.get('binning_table',  group='image_attributes', kind = 'variable')
        bin_id = binning_ids[0]

    table = f"Table_{bin_id}"
    binned_pixels = bin_data.get('bins', group=table, kind='dimension')
    binned_rows = int(binned_pixels/det_cols)
    output_data.add(name='binned_row', value=binned_rows, kind='dimension')

    return [det_rows, det_cols, binned_rows]

def check_if_binned(detector_dims):
    """
        Check if data is binned

        :param detector_dims: detector dimensions
        :type detector_dims: List

        :return is_binned: to indicated if data is binned or not
        :rtype is_binned: Boolean
    """
    is_binned = False
    det_rows = detector_dims[0]
    binned_rows = detector_dims[2]
    if det_rows != binned_rows:
        is_binned = True
    return is_binned

def add_3d_detector_image_data(data, std_data, detector_dims, output_data):
    """
        Reshape the data and add them to the netcdf output

        :param data: The data to be reshaped
        :type data: Numpy array
        :param std_data: The std data to be reshaped
        :type std data: Numpy array
        :param detector_dims: detector dimensions
        :type detector_dims: List
        :param output_data: the netcdf output
        :type output_data: dataNetcdf object
    """

    # Check if the binned row dimension is needed
    is_binned = check_if_binned(detector_dims)

    cols = detector_dims[1]
    rows = detector_dims[0]
    dimensions=('along_track','detector_row','detector_column')
    if is_binned:
        rows = detector_dims[2]
        dimensions=('along_track','binned_row','detector_column')


    data_reshaped = np.reshape(data, (data.shape[0], rows, cols))
    output_data.add(name='detector_image_3d', value=data_reshaped, group='science_data',
                    dimensions=dimensions,
                    kind='variable')
    # standard deviation data is not always present. Check
    if std_data is not None:
        std_data_reshaped = np.reshape(std_data, (std_data.shape[0], rows, cols))
        output_data.add(name='detector_stdev_3d', value=std_data_reshaped,
                        group='science_data',
                        dimensions=dimensions,
                        kind='variable')

    return


def reshape_output(output_key, config):
    """
        Reshape the output dataset
        Note: output dataset is 2D. When it is on detector level
        these dimensions are: along track and pixels.
        This second dimension is too big (det_row x det_col)
        This makes viewing difficult.
        Read in the data and reshape to 3D in case of detector level

        :param output_key: key to use to find output file in configuration
        :type output_key: String
        :param config: configuration
        :type config: Dictionary

        :return temp_output_file: file name of the output file which stores the reshaped data
        :rtype temp_output_file: String
    """
    # readin output file
    output_file = config['io'][output_key]
    output_data = dn.DataNetCDF(output_file, mode='r')

    # Get data
    data = output_data.get('detector_image',  group='science_data', kind = 'variable')
    if data is None:
        # No detector image. No need to reshape
        log.info("Not detector image present. No need to reshape and create temporary output file")
        return

    std_data = output_data.get('detector_stdev',  group='science_data', kind = 'variable')

    # Get the dimensions for the detector image and add them to the output_data
    detector_dims = detector_image_dimensions(output_data, data, config )

    add_3d_detector_image_data(data, std_data, detector_dims, output_data)

    # Remove the 2D images
    output_data.remove('detector_image', group='science_data', kind='variable')
    if output_data.find('detector_stdev', group='science_data', kind='variable') is not None:
        output_data.remove('detector_stdev', group='science_data', kind='variable')

    # write to temp output file
    temp_output_file = f"{output_file[:-3]}_temp.nc"
    output_data.write(temp_output_file)

    return temp_output_file

def reshape_output2(output_key, config):
    """
        Reshape the output dataset
        Note: output dataset is 2D. When it is on detector level
        these dimensions are: along track and pixels.
        This second dimension is too big (det_row x det_col)
        This makes viewing difficult.
        Read in the data and reshape to 3D in case of detector level

        :param output_key: key to use to find output file in configuration
        :type output_key: String
        :param config: configuration
        :type config: Dictionary

        :return temp_output_file: file name of the output file which stores the reshaped data
        :rtype temp_output_file: String
    """
    # TODO. Ik zou ook evt een check kunnen doen op bestaan van detector image of science_data ???????
    temp_output_file = None
    if is_reshape_needed(output_key, config):
        # readin output file
        output_file = config['io'][output_key]
        output_data = dn.DataNetCDF(output_file, mode='r')

        # Get data
        data = output_data.get('detector_image',  group='science_data', kind = 'variable')
        std_data = output_data.get('detector_stdev',  group='science_data', kind = 'variable')

        # Get dimensions from ckd?
        ckd_file = config['io']['ckd']
        ckd_data = dn.DataNetCDF(ckd_file, mode='r')
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
        bin_data = dn.DataNetCDF(bin_file, mode='r')
        bin_id = config['detector']['binning_table_id']
        table = f"Table_{bin_id}"
        binned_pixels = bin_data.get('bins', group=table, kind='dimension')
        binned_rows = int(binned_pixels/det_cols)
        output_data.add(name='binned_row', value=binned_rows, kind='dimension')

        # is output level with or without binning? Check original shape of data
        n_pixels = data.shape[1]
        if int(n_pixels/det_cols) == det_rows:
            # detector image
            data_reshaped = np.reshape(data, (data.shape[0], det_rows, det_cols))
            output_data.add(name='detector_image_3d', value=data_reshaped, group='science_data',
                            dimensions=('along_track','detector_row','detector_column'),
                            kind='variable')
            # standard deviation data is not always present. Check
            if std_data is not None:
                std_data_reshaped = np.reshape(std_data, (std_data.shape[0], det_rows, det_cols))
                output_data.add(name='detector_stdev_3d', value=std_data_reshaped,
                                group='science_data',
                                dimensions=('along_track','detector_row','detector_column'),
                                kind='variable')
        elif int(n_pixels/det_cols) == binned_rows:
            # binned detector image
            data_reshaped = np.reshape(data, (data.shape[0], binned_rows, det_cols))
            output_data.add(name='detector_image_3d', value=data_reshaped, group='science_data',
                            dimensions=('along_track','binned_row','detector_column'),
                            kind='variable')
            # standard deviation data is not always present. Check
            if std_data is not None:
                std_data_reshaped = np.reshape(std_data, (data.shape[0], binned_rows, det_cols))
                output_data.add(name='detector_stdev_3d', value=std_data_reshaped,
                                group='science_data',
                                dimensions=('along_track','binned_row','detector_column'),
                                kind='variable')

        # Remove the 2D images
        output_data.remove('detector_image', group='science_data', kind='variable')
        if output_data.find('detector_stdev', group='science_data', kind='variable') is not None:
            output_data.remove('detector_stdev', group='science_data', kind='variable')

        # write to temp output file
        temp_output_file = f"{output_file[:-3]}_temp.nc"
        output_data.write(temp_output_file)

    return temp_output_file

def get_file_name(config, step):
    """
        Determine output/input file path.
        When it is a nominal run it is just base_path
        When it concerns a scenarion and a step is required to be rerun
        for this scenario it is base_path/scenarios/scenario_subdir

        :param: config: the configuration
        :type:  config: Dictionary
        :param: step: the step of the E2E which needs to be run
        :type:  step: String

        :return file_path: Path to requested file
        :rtype: file_path: String
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

def get_specific_config(orig_config, kind):
    """
        Obtain module specific configuration from full configuration
        Possible module names: gm, sgm, im, l1al1b, l1l2, pam

        :param: orig_config: original full configuration
        :type: orig_config: Dictionary
        :param: kind: the name of the module
        :type: kind: String

        :return specific_config: Configuration specific for given module
        :rtype specific_config: Dictionary
    """
    if kind not in orig_config:
        error_msg = f"Unknown module {kind}. Unable to fetch corresponding config. Exiting!"
        log.error(error_msg)
        sys.exit(error_msg)

    # Get the module specific settings
    specific_config = orig_config[kind]
    # Get the header
    specific_config['header'] = orig_config['header']
    # Get the module IO settings
    specific_config['io'] = {}

    if kind == 'gm':
        # Combine path and file name
        output_path = get_file_name(orig_config, 'gm')
        specific_config['io']['gm'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'sgm':
        # Combine path and file name
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm_rad'] = os.path.join(output_path, orig_config['io']['sgm_rad'])
        specific_config['io']['sgm_atm_raw'] = os.path.join(output_path, \
            orig_config['io']['sgm_atm_raw'])
        specific_config['io']['sgm_atm'] = os.path.join(output_path, orig_config['io']['sgm_atm'])
        specific_config['io']['gm'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'im':
        specific_config['io']['binning_table'] = orig_config['io']['binning_table']
        specific_config['io']['ckd'] = orig_config['io']['ckd_im']

        do_python = specific_config['do_python']
        # Output of IM
        output_path = get_file_name(orig_config, 'im')
        l1a_file_name = orig_config['io']['l1a']
        if do_python:
            l1a_file_name_python = f"{l1a_file_name[0:-3]}_python.nc"
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name_python)
            specific_config['io']['im_algo_output'] = os.path.join(output_path, \
                orig_config['io']['im_algo_output'])
        else:
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name)

        # Input to IM
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm'] = os.path.join(output_path, orig_config['io']['sgm_rad'])

        # To add info from geometry file
        output_path = get_file_name(orig_config, 'gm')
        specific_config['io']['geometry'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'l1al1b':

        # Output to L1B
        output_path = get_file_name(orig_config, 'l1al1b')
        l1b_file_name = orig_config['io']['l1b']

        do_python = specific_config['do_python']
        if do_python:
            l1b_file_name_python = f"{l1b_file_name[0:-3]}_python.nc"
            specific_config['io']['l1b'] = os.path.join(output_path, l1b_file_name_python)
            specific_config['io']['l1b_algo_output'] = os.path.join(output_path, \
                orig_config['io']['l1b_algo_output'])
        else:
            specific_config['io']['l1b'] = os.path.join(output_path, l1b_file_name)

        # Input to L1B
        output_path = get_file_name(orig_config, 'im')
        l1a_file_name = orig_config['io']['l1a']

        # Note: if IM was run using the python code, the l1a python output needs to be readin
        do_python = orig_config['im']['do_python']
        if do_python:
            l1a_file_name_python = f"{l1a_file_name[0:-3]}_python.nc"
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name_python)
        else:
            specific_config['io']['l1a'] = os.path.join(output_path, l1a_file_name)

        specific_config['io']['binning_table'] = orig_config['io']['binning_table']
        specific_config['io']['ckd'] = orig_config['io']['ckd']

        output_path = get_file_name(orig_config, 'gm')
        specific_config['io']['geometry'] = os.path.join(output_path, orig_config['io']['gm'])

    elif kind == 'l1l2':
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
        specific_config['isrf'] = orig_config['im']['isrf']

    elif kind == 'pam':
        output_path = get_file_name(orig_config, 'sgm')
        specific_config['io']['sgm_atm'] = os.path.join(output_path, orig_config['io']['sgm_atm'])

        output_path = get_file_name(orig_config, 'l1l2')
        specific_config['io']['l2'] = os.path.join(output_path, orig_config['io']['l2'])

    else:
        error_msg = f"Unknown module kind {kind}. Unable to fetch corresponding IO config. Exiting!"
        log.error(error_msg)
        sys.exit(error_msg)

    return specific_config

def add_module_specific_attributes(config, attribute_dict, step):
    """
        Add the module specific settings to the attribute_dict

        :param: config: configuration
        :type:  config: Dictionary
        :param: attribute_dict: dictionary of file attributes
        :type:  attribute_dict: Dictionary
        :param: step: name of the module
        :type:  step: String

        :return: new_attribute_dict: new dictionary of file attributes
        :rtype:  new_attribute_dict: Dictionary
    """
    new_attribute_dict = attribute_dict.copy()

    for key, value in config.items():
        new_key = f"{step}_{key}"
        new_attribute_dict[new_key] = value

    return new_attribute_dict


def build(config, step, cfg_path, attribute_dict):
    """
        Run E2E processor.
        :param: config: configuration file containing the settings for the different
                steps in the E2E processor
        :type:  config: Dictionary
        :param: step: indicating which step in the E2E processor to run.
        :type:  step: String
        :param: cfg_path: configuration path
        :type:  cfg_path: String
        :param: attribute_dict: Dictionary with attributes to be added to main of output netCDF file
        :type:  attribute_dict: Dictionary
    """

    # Make a copy of the full config file.
    configuration = config.copy()

    if step in ['gm','all']:

        gm_config = get_specific_config(configuration, 'gm')
        attribute_dict = add_module_specific_attributes(gm_config, attribute_dict, 'gm')
        e2e_module = importlib.import_module("gm.gm")
        e2e_module.geometry_module(gm_config)
        # add attributes to the output file
        Utils.add_attributes_to_output(gm_config['io']['gm'], attribute_dict)

    if step in ['sgm','all']:
        sgm_config = get_specific_config(configuration, 'sgm')
        attribute_dict = add_module_specific_attributes(sgm_config, attribute_dict, 'sgm')
        e2e_module = importlib.import_module("sgm.sgm_no2")
        e2e_module.scene_generation_module_nitro(sgm_config)
        Utils.add_attributes_to_output(sgm_config['io']['sgm_rad'], attribute_dict)
        Utils.add_attributes_to_output(sgm_config['io']['sgm_atm'], attribute_dict)
        Utils.add_attributes_to_output(sgm_config['io']['sgm_atm_raw'], attribute_dict)

    if step in ['im','all']:
        im_config = get_specific_config(configuration, 'im')
        attribute_dict = add_module_specific_attributes(im_config, attribute_dict, 'im')
        # write config to temp yaml file with IM values filled in
        im_config_file = f"{cfg_path}/im_config_temp.yaml"
        with open(im_config_file,'w') as outfile:
            yaml.dump(im_config, outfile)

        if not im_config['do_python']:
            # Run C++ IM code
            # Need to call C++ using the IM specific config_file
            # TODO: Should check be set to True
            subprocess.run(["im/build/tango_im.x", im_config_file], check=False)

# With updates to l1_measurement.cpp there is no need for reshaping data.
#            # output dataset is 2D. In case of detector image (in case of some inbetween steps
#            # and of the final output) the second dimension is detector_pixels which is too
#            # large to view.
#            # Need to reshape to 3D to be able to make sense of this.
#            temp_output_file = reshape_output('l1a', im_config)

            # Add attributes to output file
#            if temp_output_file is not None:
#                Utils.add_attributes_to_output(temp_output_file, attribute_dict)
            Utils.add_attributes_to_output(im_config['io']['l1a'], attribute_dict)
        else:
            # run Python code
            e2e_module = importlib.import_module("im.Python.instrument_model")
            e2e_module.instrument_model(im_config, attribute_dict)
            # No need to reshape because Python output is already 3D

    if step in ['l1al1b','all']:

        l1b_config = get_specific_config(configuration, 'l1al1b')
        attribute_dict = add_module_specific_attributes(l1b_config, attribute_dict, 'l1al1b')
        # write config to temp yaml file with L1B values filled in
        l1b_config_file = f"{cfg_path}/l1b_config_temp.yaml"
        with open(l1b_config_file,'w') as outfile:
            yaml.dump(l1b_config, outfile)

        if not l1b_config['do_python']:

            # Need to call C++ with L1B specific config file
            subprocess.run(["l1al1b/build/tango_l1b.x", l1b_config_file], check=False)

# With updates to l1_measurement.cpp there is no need for reshaping data.
#            # output dataset is 2D. In case of detector image (in case of inbetween step)
#            # the second dimension is detector_pixels which is too large to view.
#            # Need to reshape to 3D to be able to make sense of this.
#            temp_output_file = reshape_output('l1b', l1b_config)

            # Add attributes to output file
#            if temp_output_file is not None:
#                Utils.add_attributes_to_output(temp_output_file, attribute_dict)
            Utils.add_attributes_to_output(l1b_config['io']['l1b'], attribute_dict)
        else:
            # run Python code
            e2e_module = importlib.import_module("l1al1b.Python.l1b_processor")
            e2e_module.l1b_processor(l1b_config, attribute_dict)

    if step in ['l1l2','all']:
        l2_config = get_specific_config(configuration, 'l1l2')
        attribute_dict = add_module_specific_attributes(l2_config, attribute_dict, 'l1l2')
        e2e_module = importlib.import_module("l1l2.l1bl2_no2")
        e2e_module.l1bl2_no2(l2_config)
        Utils.add_attributes_to_output(l2_config['io']['l2'], attribute_dict)

    if step in ['pam','all']:
        pam_config = get_specific_config(configuration, 'pam')
        attribute_dict = add_module_specific_attributes(pam_config, attribute_dict, 'pam')
        e2e_module = importlib.import_module("pam.pam")
        e2e_module.pam_nitro(pam_config)

if __name__ == "__main__":

    startTime = time.time()

    # Get input arguments
    e2e_cfg_file, e2e_step =  cmdline(sys.argv[1:])
    e2e_cfg_path, filename = os.path.split(e2e_cfg_file)

    # Get configuration info
    e2e_config = Utils.get_config(e2e_cfg_file)
    e2e_config['header']['path'] = e2e_cfg_path

    if e2e_config['log_level'] == 'debug':
        log.setLevel(_logging.DEBUG)
    else:
        log.setLevel(_logging.INFO)

    if not e2e_config['show_progress']:
        os.environ['TQDM_DISABLE']='1'


    # Get information (like git hash and config file name and version (if available)
    # that will be added to the output files as attributes
    main_attribute_dict = Utils.get_main_attributes(e2e_config)

    build(e2e_config, e2e_step, e2e_cfg_path, main_attribute_dict)

    log.info(f'E2E calculation finished in {timedelta(seconds=time.time()-startTime)}')
