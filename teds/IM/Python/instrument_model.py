import logging
import sys
import time
import numpy as np

import teds.lib.lib_utils as Utils
import teds.lib.data_netcdf.data_netcdf as dn

from teds.IM.Python.datasets import Datasets
from teds.IM.Python.algos.algo_wavemap import Wavemap

from teds.IM.Python.input.input_yaml import Input_yaml
from teds.IM.Python.input.input_netcdf import Input_netcdf

from teds.IM.Python.algos.algo_isrf import ISRF
from teds.IM.Python.algos.algo_draw_on_detector import Draw_On_Detector
from teds.IM.Python.algos.algo_simple_regrid import Simple_Regrid
from teds.IM.Python.algos.algo_prnu import PRNU
from teds.IM.Python.algos.algo_dark_current import Dark_Current


def get_input_data(logger, config):
    """
        Get the input data.
        ckd, proctable, radiance, wavelength,...
        return them in input_data dictionary
    """

    input_data = Datasets(logger, 'input_data')

    # Proctable
    proctable_file = config['proctable']['file']
    proctable_input = Input_yaml(logger, proctable_file)
    proctable = proctable_input.read()
    input_data.add_container('proctable', proctable)

    # CKD
    ckd_file = config['io']['ckd_im']
    ckd_input = Input_netcdf(logger,ckd_file)
    ckd = ckd_input.read()
    #Note: ckd is an data_netcdf object
    # It contains several groups and datasets
    input_data.add_container('ckd', ckd)

    input_file = config['sgm_rad_file']
    data_input = Input_netcdf(logger,input_file)
    #Note: data_input is an data_netcdf object
    # It contains several groups and datasets
    input_datasets = data_input.read()
    input_data.add_container('measurement', input_datasets)

    # Find the wavelength map corresponding to this temperature.
    # Not yet implemented
    wavemapping = Wavemap(logger)
    wavemapping.check_input(input_data)
    wavemapping.execute(input_data)
    wavemap = wavemapping.get_data()
    input_data.add_container('extra', {'wavemap':wavemap})

    # Also add configuration settings to input_data
    input_data.add_container('config', config)

    return input_data

def  add_main_attributes(out_data, attributes):
    """
        Add attributes to output file
    """
    for name, value in attributes.items():
        out_data.add(name, value=value, kind='attribute')
    return

def initialize_output(input_data, algo_list, dimensions, main_attributes):
    """
        Initialize the output files.
        Add dimensions and output dataset
    """
    output_datasets = Datasets(logger, 'output_data')

    # Create and initialize the final output
    output_file_name = input_data.get_dataset('l1a', c_name='config', group='io')
    output = dn.DataNetCDF(logger, output_file_name)
    add_main_attributes(output, main_attributes)

    # Figure out which dimensions should be in the final output file
    # If algo rithm Draw_On_Detector is in the algo_list, output id on detector pixels
    # If ISRF in algo list it is mixed
    # Is this check sufficient?
    # Check if draw_on_detector algo is in algo list:
    output.add(name='scanline', value=dim_alt, kind='dimension') 
    if 'Draw_On_Detector' in algo_list:
        # we have a detector image so detector dimensions can be used in output
        output.add(name='row', value=dim_spat, kind='dimension') 
        output.add(name='col', value=dim_spec, kind='dimension') 
        output_data = np.zeros((dim_alt, dim_spat, dim_spec))
        output.add(name='measurement', value=output_data, dimensions=('scanline','row','col'), kind='variable')
    elif 'ISRF' in algo_list:
        print(f"Creating ISRF dimensions for final output")
        # detector columns combined with actrack dimension
        output.add(name='act', value=dim_act, kind='dimension') 
        output.add(name='col', value=dim_spec, kind='dimension') 
        output_data = np.zeros((dim_alt, dim_act, dim_spec))
        output.add(name='measurement', value=output_data, dimensions=('scanline','act','col'), kind='variable')
    else:
        # ISRF AND Draw_On_Detector not in algo_list. Any other algos are on detector pixels. Can not be applied.
        # Can not continue
        error_message = " Both the ISRF and the Draw_On_Detector algos are not involved. All other algos run on detector pixels and can not be applied. Can not continue"
        logger.error(error_message)
        sys.exit(error_message)

    output_datasets.add_container('final',output)

    # Now add output for different steps
    for algo_name in algo_list:

        # These are the inbetween output files
        # Maybe add some switch if we want then or not.
        algo_output = input_data.get_dataset('im_algo_output', c_name='config', group='io')
        algo_file = algo_output.format(algo_name=algo_name)
        output_algo = dn.DataNetCDF(logger, algo_file)
        add_main_attributes(output_algo, main_attributes)

        # Add dimensions to the inbetween output files
        # In principle output dimensions for ISRF is mixed
        # All others should be detector pixels 
        output_algo.add(name='scanline', value=dimensions['dim_alt'], kind='dimension') 
        if algo_name == 'ISRF':
            # output not yet on full detector dimensions
            output_algo.add(name='act', value=dimensions['dim_act'], kind='dimension') 
            output_algo.add(name='col', value=dimensions['dim_spec'], kind='dimension') 
            # Create output dataset with zeros. For now named it measurement. Check what name of output dataset should be
            # When running over images this dataset will be updated
            output_algo_data = np.zeros((dimensions['dim_alt'], dimensions['dim_act'], dimensions['dim_spec']))
            output_algo.add(name='measurement', value=output_algo_data, dimensions=('scanline','act','col'), kind='variable')
        else:
            output_algo.add(name='row', value=dimensions['dim_spat'], kind='dimension') 
            output_algo.add(name='col', value=dimensions['dim_spec'], kind='dimension') 
            # Create output dataset with zeros. For now named it measurement. Check what name of output dataset should be
            # When running over images this dataset will be updated
            output_algo_data = np.zeros((dimensions['dim_alt'], dimensions['dim_spat'], dimensions['dim_spec']))
            output_algo.add(name='measurement', value=output_algo_data, dimensions=('scanline','row','col'), kind='variable')

        # Maybe add group original to store the input data????

        output_datasets.add_container(algo_name,output_algo)

    return output_datasets

def instrument_model(config, logger, main_attributes):
    """
        Instrument model.
        Apply effects to create L0 data
    """

    start_time_IM = time.perf_counter()

    # Get the input data into list of data containers
    input_data = get_input_data(logger, config)

    # Get the list of algorithms from proctable file
    scenario = input_data.get_dataset('scenario', c_name='config', group='proctable')
    algo_list = input_data.get_dataset(scenario, c_name='proctable')
    logger.info(f"NOW algo_list: {algo_list}")

    # Get the input data. Note: possibly this needs to be updated when SGM is updated.
    radiance_data = input_data.get_dataset('radiance', c_name='measurement', kind='variable')

    #Note: this is the along track dimension
    n_images = radiance_data.shape[0]

    dim_alt = n_images
    dim_spat = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
    dim_spec = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
    dim_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')
    logger.info(f"Found dim_spec: {dim_spec} and dim_spat: {dim_spat} and dim_act: {dim_act}")
    dimensions = {'dim_act':dim_act,'dim_spec': dim_spec, 'dim_spat': dim_spat. 'dim_alt': dim_alt}

    # Get the output data into list of data containers
    output_datasets = initialize_output(input_data, algo_list, dimensions, main_attributes)

    # Loop over images
    for img in range(n_images):

        start_time_image = time.perf_counter()

        image = radiance_data[img,:,:]
        # Have to temporary store the image somewhere for acces later on in algos.
        input_data.add_container('work', {'image': image})

        for algo_name in algo_list:
            start_time_algo = time.perf_counter()

            # TODO in principle this works. Check if more complicated stuff also works
            logger.debug(f"IMG: {img} running algo {algo_name}")

            algo = globals()[algo_name](logger)
            algo.check_input(input_data)

            logger.debug(f"BEFORE RUNNING ALGO SHAPE IMAGE: {image.shape}")

            algo.execute( input_data)
            image = algo.get_data()

            logger.debug(f"AFTER RUNNING ALGO SHAPE IMAGE: {image.shape}")

            logger.debug(f"Updating measurement dataset for container {algo_name} with image nr {img}")
            output_datasets.update_dataset('measurement', data=image, c_name=algo_name, img=img)

            input_data.update_dataset('image', 'work',image)
            this_dataset = output_datasets.get_dataset('measurement', c_name=algo_name, kind='variable')

            logger.debug(f"{algo_name} SHAPE : {this_dataset.shape}")

            end_time_algo = time.perf_counter()
            logger.info(f"Processing algorithm {algo_name} for image {img} took {(end_time_algo - start_time_algo):.6f}s")

        end_time_image = time.perf_counter()
        logger.info(f"Processing image {img} took {(end_time_image - start_time_image):.6f}seconds")

    # Get dataset from last algo
    final_dataset = output_datasets.get_dataset('measurement', c_name=algo_list[-1], kind='variable')
    # update dataset from final container
    output_datasets.update_dataset('measurement', data=this_dataset, c_name='final')

    output_datasets.write()
    end_time_IM = time.perf_counter()
    logger.info(f"Instrument Model Processing took {(end_time_IM - start_time_IM):.6f}seconds")

    logger.info("Succes")

if __name__ == '__main__' :

    # Get logger for IM
    im_logger = Utils.get_logger()

    # Get configuration info
    cfgFile = sys.argv[1]
    config = Utils.getConfig(im_logger, cfgFile)

    # Get information (like git hash and config file name and version (if available) 
    # that will be added to the output file as attributes
    main_attribute_dict = Utils.get_main_attributes(config, config_attributes_name='IM_configuration')
    # Also get some other attributes?

    # Get started
    instrument_model(config, im_logger, main_attribute_dict)


