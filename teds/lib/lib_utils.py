import os, sys
import logging
import yaml
import subprocess
import lib.data_netcdf.data_netcdf as dn

def get_logger():
    """
       Gets or creates a logger
    """

    log_level = logging.INFO
    log_format = '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s'
    log_formatter = logging.Formatter(log_format)
    date_format = '%d/%m/%Y %H:%M:%S'

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    c_handler = logging.StreamHandler(sys.stdout)
    c_handler.setFormatter(log_formatter)
    logger.addHandler(c_handler)

    return logger

def getConfig(logger, cfgFile):
    """
        Get the config information from the configuration file
       - logger: Reference to the program logger
       - cfgFile: configuration file
       return:
       - configuration: configuration info
    """
    cfg_path, filename = os.path.split(cfgFile)
    stream =  open(cfgFile, 'r')
    config = yaml.safe_load(stream)
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

    logger.info(config_string)

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
        config.pop('header')
#    attribute_dict['E2E_configuration'] = str(config)
    attribute_dict[config_attributes_name] = str(config)
    #TODO: Add other information that might be handy to have in the attributes of the netCDF output file.

    return attribute_dict

def add_attributes_to_output(logger, output_file, attribute_dict):
    """
        Add attributes to the output file
    """

    out_data = dn.DataNetCDF(logger, output_file, mode='r')
    for name, value in attribute_dict.items():
        out_data.add(name, value=value, kind='attribute')
    out_data.write()
    return

