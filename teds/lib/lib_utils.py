import os, sys
import logging
import subprocess
import netCDF4 as nc

from teds.im.Python.input.input_yaml import Input_yaml
# import teds.lib.data_netcdf.data_netcdf as dn
from teds import log

def get_logger():
    """
       Initialise logger with name 'E2E'.
       Logger can be obtained in other modules with:
            logger = logging.getLogger('E2E')
    """

    log_level = logging.INFO
    log_format = '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s'
    log_formatter = logging.Formatter(log_format)

    logger = logging.getLogger('E2E')
    logger.setLevel(log_level)

    c_handler = logging.StreamHandler(sys.stdout)
    c_handler.setFormatter(log_formatter)
    logger.addHandler(c_handler)

    return logger

def getConfig(cfgFile):
    """
        Get the config information from the configuration file
       - logger: Reference to the program logger
       - cfgFile: configuration file
       return:
       - configuration: configuration info
    """
    cfg_path, filename = os.path.split(cfgFile)

    config_input = Input_yaml(log, cfgFile)
    config = config_input.read()
    # print(f"{config_input}")

#    stream =  open(cfgFile, 'r')
#    config = yaml.safe_load(stream)
    if 'header' in config.keys():
        config['header']['path']=cfg_path
    else:
        config['header'] = {'path':cfg_path, 'file_name': filename, 'version': 'not set'}

    #TODO do we need the below capability?
    # Fill in variables in the configuration
    for key in config:
        if isinstance(config[key], str):
            config[key] = config[key].format(**config)
        if isinstance(config[key], list):
            for i, item in enumerate(config[key]):
                if isinstance(item, str):
                    config[key][i] = config[key][i].format(**config)

    config_string = "Configuration used in this step of the analysis:\n"
    for key in config:
        config_string += f"    {key} = {config[key]}\n"


    config_string = config_input.print()
    log.info(config_string)

    return config

def get_main_attributes(config, config_attributes_name='E2E_configuration'):
    """
        Define some info to add as attributes to the main of the ouput netCDF file
    """

    attribute_dict = {}

    attribute_dict['git_hash'] = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip())
    attribute_dict['git_hash_short'] = str( subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip())

    # If header information exists also put it as attribute in output file
    if 'header' in config.keys():
        attribute_dict['cfg_file'] = config['header']['file_name']
        attribute_dict['cfg_version'] = config['header']['version']
        attribute_dict['cfg_path'] = config['header']['path']
#        config.pop('header')
    attribute_dict[config_attributes_name] = str(config)
    if 'scenario' in config:
        scenario_settings = config['scenario']
        attribute_dict['scenario_title'] = scenario_settings['title']
        attribute_dict['scenario_description'] = scenario_settings['description']
        attribute_dict['scenario_subdir'] = scenario_settings['subdir']
        attribute_dict['scenario_steps'] = scenario_settings['steps']

    return attribute_dict

def add_attributes_to_output(output_file, attribute_dict):
    """
        Add attributes to the output file
    """

    # out_data = dn.DataNetCDF(logger, output_file, mode='r')
    # for name, value in attribute_dict.items():
    #     out_data.add(name, value=str(value), kind='attribute')
    # out_data.write()
    # return

    with nc.Dataset(output_file, 'a') as file:
        for name, value in attribute_dict.items():
            # TODO The line below does not seem to be working
            setattr(file, name, str(value))
    return
