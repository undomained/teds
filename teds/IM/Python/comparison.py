import logging
import sys
import numpy as np

import teds.lib.lib_utils as Utils
import teds.lib.data_netcdf.data_netcdf as dn
from teds.IM.Python.datasets import Datasets
from teds.IM.Python.input.input_yaml import Input_yaml
from teds.IM.Python.input.input_netcdf import Input_netcdf

def find_dataset(logger, input_data, group_name, dataset_name, data_file_1, data_file_2):

    var_1 = input_data.get_dataset(dataset_name, c_name='data_1', group=group_name, kind='variable')
    var_2 = input_data.get_dataset(dataset_name, c_name='data_2', group=group_name, kind='variable')

    print(f"GROUP: {group_name} and dataset: {dataset_name}, var_1: {var_1} and var_2: {var_2}")

    if var_1 is None and var_2 is not None:
        # dataset name could be detector_image_3d in case of C++ file
        ds_updated = f'{dataset_name}_3d'
        var_1 = input_data.get_dataset(ds_updated, c_name='data_1', group=group_name, kind='variable')
        if var_1 is None:
            # Some other problem, Bailing out
            error_msg = f"Variable {dataset_name} in group {group_name} found in file ({data_file_2}), but not found in file ({data_file_1}). No comparison possible!"
            logger.error(error_msg)
            sys.exit(error_msg)
    if var_2 is None and var_1 is not None:
        ds_updated = '{dataset_name}_3d'
        var_2 = input_data.get_dataset(ds_updated, c_name='data_2', group=group_name, kind='variable')
        if var_2 is None:
            # Some other problem, Bailing out
            error_msg = f"Variable {dataset_name} in group {group_name} found in file ({data_file_1}), but not found in file ({data_file_1}). No comparison possible!"
            logger.error(error_msg)
            sys.exit(error_msg)

    return var_1, var_2

def compare_data(logger, data_file_1, data_file_2):

    data_1_input = Input_netcdf(logger,data_file_1)
    data_1 = data_1_input.read()
    print(f"DATA_1: {data_1}")
    data_2_input = Input_netcdf(logger,data_file_2)
    data_2 = data_2_input.read()
    print(f"DATA_2: {data_2}")

    input_data = Datasets(logger, 'input_data')
    input_data.add_container('data_1', data_1)
    input_data.add_container('data_2', data_2)
    print(f"Added data_1 and data_2 to input_data")

    group_name = 'science_data'
    dataset_name = 'detector_image'
    logger.info(f"Trying to find {dataset_name} in group {group_name}")
    var_1, var_2 = find_dataset(logger, input_data, group_name, dataset_name, data_file_1, data_file_2)

    if var_1 is None and var_2 is None:
        # try observation data
        group_name = 'observation_data'
        dataset_name = "i"
        logger.info(f"Trying to find {dataset_name} in group {group_name}")
        var_1, var_2 = find_dataset(logger, input_data, group_name, dataset_name, data_file_1, data_file_2)

        if var_1 is None and var_2 is None:
            # Both files do not contain detector image and no observation data. Nothing to compare. 
            error_msg = f"Both files do not contain detector image and no observation data. Nothing to compare."
            logger.error(error_msg)
            return

    print(f"Shape: var_1: {var_1.shape}, var_2: {var_2.shape}")
    # If we get here we have data to compare
    difference = np.subtract(var_1,var_2)
    division = np.divide(var_1,var_2)

    output_file_name = 'comparison.nc'
    nc_output = dn.DataNetCDF(logger, output_file_name)

    dims = var_1.shape
    if group_name == 'science_data':
        nc_output.add(name='detector_image', value=dims[0], kind='dimension')
        nc_output.add(name='row', value=dims[1], kind='dimension')
        nc_output.add(name='col', value=dims[2], kind='dimension')
        nc_output.add(name='science_data', kind='group')
        nc_output.add(name='original_data_1', value=var_1, dimensions=('detector_image','row','col'), group='science_data', kind='variable')
        nc_output.add(name='original_data_2', value=var_2, dimensions=('detector_image','row','col'), group='science_data', kind='variable')
        nc_output.add(name='difference', value=difference, dimensions=('detector_image','row','col'), group='science_data', kind='variable')
        nc_output.add(name='division', value=division, dimensions=('detector_image','row','col'), group='science_data', kind='variable')
    elif group_name == 'observation_data':
        nc_output.add(name='along_track', value=dims[0], kind='dimension')
        nc_output.add(name='across_track', value=dims[1], kind='dimension')
        nc_output.add(name='wavelength', value=dims[2], kind='dimension')
        nc_output.add(name='observation_data', kind='group')
        nc_output.add(name='original_data_1', value=var_1, dimensions=('along_track','across_track','wavelength'), group='observation_data', kind='variable')
        nc_output.add(name='original_data_2', value=var_2, dimensions=('along_track','across_track','wavelength'), group='observation_data', kind='variable')
        nc_output.add(name='difference', value=difference, dimensions=('along_track','across_track','wavelength'), group='observation_data', kind='variable')
        nc_output.add(name='division', value=division, dimensions=('along_track','across_track','wavelength'), group='observation_data', kind='variable')

    nc_output.write()

    logger.info("Succes")

if __name__ == '__main__' :

    comp_logger = Utils.get_logger()

    data_1, data_2 = sys.argv[1:]

    comp_logger.info(f"Comparing files {data_1} and {data_2}")
    compare_data(comp_logger, data_1, data_2)
