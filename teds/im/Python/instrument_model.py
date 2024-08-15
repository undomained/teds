# pylint: disable=unused-import, fixme, logging-fstring-interpolation
"""Providing code to run the Tango instrument module"""

import sys
import time
import numpy as np

import teds.lib.lib_utils as Utils
import teds.lib.data_netcdf.data_netcdf as dn

from teds import log
from teds.run_E2E import get_specific_config, add_module_specific_attributes

from teds.im.Python.datasets import Datasets
from teds.im.Python.algos.algo_wavemap import Wavemap

from teds.im.Python.input.input_yaml import Input_yaml
from teds.im.Python.input.input_netcdf import Input_netcdf

from teds.im.Python.algos.algo_isrf import ISRF
from teds.im.Python.algos.algo_radiometric import Radiometric
from teds.im.Python.algos.algo_draw_on_detector import Draw_On_Detector
from teds.im.Python.algos.algo_simple_regrid import Simple_Regrid
from teds.im.Python.algos.algo_prnu import PRNU
from teds.im.Python.algos.algo_dark_current import Dark_Current
from teds.im.Python.algos.algo_noise import Noise
from teds.im.Python.algos.algo_dark_offset import Dark_Offset
from teds.im.Python.algos.algo_coadding import Coadding
from teds.im.Python.algos.algo_binning import Binning
from teds.im.Python.algos.algo_adc import ADC

def get_im_config(config):
    """
        Get IM specific settings and move them one level up in config
    """
    im_settings = config['im']
    for key, value in im_settings.items():
        config[key] = value
    return config

def get_input_data(config):
    """
        Get the input data.
        ckd, proctable, radiance, wavelength,...
        return them in input_data dictionary
    """

    input_data = Datasets(log, 'input_data')

    # Proctable
    proctable_file = config['proctable']['file']
    proctable_input = Input_yaml(log, proctable_file)
    proctable = proctable_input.read()
    input_data.add_container('proctable', proctable)

    # CKD
    ckd_file = config['io']['ckd']
    ckd_input = Input_netcdf(log,ckd_file)
    ckd = ckd_input.read()
    #Note: ckd is an data_netcdf object
    # It contains several groups and datasets
    input_data.add_container('ckd', ckd)

    # Binning Tables
    binning_file = config['io']['binning_table']
    print(f"binning file: {binning_file}")
    binning_tables = Input_netcdf(log,binning_file)
    binning_table_datasets = binning_tables.read()
    input_data.add_container('binning', binning_table_datasets)

    input_file = config['io']['sgm']
    # input_file = config['io']['l1b']
    data_input = Input_netcdf(log,input_file)
    #Note: data_input is an data_netcdf object
    # It contains several groups and datasets
    input_datasets = data_input.read()
    input_data.add_container('measurement', input_datasets)

    # Find the wavelength map corresponding to this temperature.
    # Not yet implemented
    wavemapping = Wavemap(log)
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
        out_data.add(name=name, value=str(value), kind='attribute')
    return

def is_detector_image(algo_name):
    """
        Determine if th output are detector images or observation data
    """
    is_a_detector_image = False
    # TODO. This is now hard coded. Need to do this differently
    detector_image_algos = ['Draw_On_Detector', 'PRNU','Non_Linearity','Straylight',
                            'Dark_Current', 'Noise', 'Dark_Offset','Coadding','Binning','ADC']
    if algo_name in detector_image_algos:
        is_a_detector_image = True
    return is_a_detector_image

def create_detector_image_output(nc_output, dimensions, image_attribute_data,
    is_binned=False, in_between=False
):
    """
        Create detector image dimensions and dataset and attributes
    """
    # Detector image
    nc_output.add(name='detector_image', value=dimensions['dim_alt'], kind='dimension')
    nc_output.add(name='col', value=dimensions['dim_spec'], kind='dimension')
    # Binned rows?
    if is_binned:
        # output is binned
        nc_output.add(name='row', value=dimensions['dim_binned_rows'], kind='dimension')
        output_data = np.zeros((dimensions['dim_alt'], dimensions['dim_binned_rows'],
                                dimensions['dim_spec']))
    else:
        nc_output.add(name='row', value=dimensions['dim_spat'], kind='dimension')
        output_data = np.zeros((dimensions['dim_alt'], dimensions['dim_spat'],
                                dimensions['dim_spec']))

    nc_output.add(name='science_data', kind='group')
    nc_output.add(name='detector_image', value=output_data,
                  dimensions=('detector_image','row','col'),
                  group='science_data', kind='variable')

    nc_output.add(name='name', value='detector_images', var='detector_image',
                  group='science_data', kind='attribute')
    nc_output.add(name='units', value='counts', var='detector_image',
                  group='science_data',kind='attribute')

    if in_between:
        nc_output.add(name='measurement', value=output_data,
                      dimensions=('detector_image','row','col'), kind='variable')

    nc_output.add(name='image_attributes', kind='group')
    nc_output.add(name='binning_table', value=image_attribute_data['binning_data'],
                  dimensions=('detector_image'), group='image_attributes', kind='variable')
    nc_output.add(name='exposure_time', value=image_attribute_data['expt_data'],
                  dimensions=('detector_image'), group='image_attributes', kind='variable')
    nc_output.add(name='nr_coadditions', value=image_attribute_data['coadd_data'],
                  dimensions=('detector_image'), group='image_attributes', kind='variable')
    nc_output.add(name='image_time', value=image_attribute_data['image_time_data'],
                  dimensions=('detector_image'), group='image_attributes', kind='variable')

    # TODO need to add more attributes?
    # TODO also add standard deviation data (detector_stdev)?

    return

def create_observation_data_output(nc_output, dimensions, image_attribute_data, in_between=False):
    """
        Create observation data dimensions and dataset and attributes
    """
    # Observation data
    nc_output.add(name='along_track', value=dimensions['dim_alt'], kind='dimension')
    nc_output.add(name='across_track', value=dimensions['dim_act'], kind='dimension')
    nc_output.add(name='wavelength', value=dimensions['dim_spec'], kind='dimension')
    output_data = np.zeros((dimensions['dim_alt'], dimensions['dim_act'], dimensions['dim_spec']))

    nc_output.add(name='observation_data', kind='group')
    nc_output.add(name='i', value=output_data,
                  dimensions=('along_track','across_track','wavelength'),
                  group='observation_data', kind='variable')
    nc_output.add(name='name', value='radiance', var='i',
                  group='observation_data', kind='attribute')
    nc_output.add(name='units', value='ph nm-1 s-1 sr-1 m-2', var='i',
                  group='observation_data', kind='attribute')

    if in_between:
        nc_output.add(name='measurement', value=output_data,
                      dimensions=('along_track','across_track','wavelength'), kind='variable')

    nc_output.add(name='image_attributes', kind='group')
    nc_output.add(name='binning_table', value=image_attribute_data['binning_data'],
                  dimensions=('along_track'), group='image_attributes', kind='variable')
    nc_output.add(name='exposure_time', value=image_attribute_data['expt_data'],
                  dimensions=('along_track'), group='image_attributes', kind='variable')
    nc_output.add(name='nr_coadditions', value=image_attribute_data['coadd_data'],
                  dimensions=('along_track'), group='image_attributes', kind='variable')
    nc_output.add(name='image_time', value=image_attribute_data['image_time_data'],
                  dimensions=('along_track'), group='image_attributes', kind='variable')

    # TODO add other attributes?????? or need to add standard deviation data (i_stdev)?
    # Or need to add group sensor_bands and wavelength data?
    return


def initialize_output(kind, input_data, algo_list, dimensions, attributes):
    """
        Initialize the final output file.
        Add dimensions and output dataset
        kind is l1a or l1b
    """

    output_datasets = Datasets(log, 'output_data')

# Final output file
    # Create and initialize the final output
    output_file_name = input_data.get_dataset(kind, c_name='config', group='io')
    output = dn.DataNetCDF(log, output_file_name)
    add_main_attributes(output, attributes)

    # Output file also needs image_time, binning_table (=id) nr_coadditions, exposure_time
    # all as fct of image_nr
    n_images = dimensions['dim_alt']
    nr_coads = input_data.get_dataset('nr_coadditions', c_name='config', group='detector')
    expt = input_data.get_dataset('exposure_time', c_name='config', group='detector')
    binning = input_data.get_dataset('binning_table_id', c_name='config', group='detector')
    image_attribute_data = {}
    image_attribute_data['binning_data'] = binning* np.ones((n_images,))
    image_attribute_data['expt_data'] = expt* np.ones((n_images,))
    image_attribute_data['coadd_data'] = nr_coads* np.ones((n_images,))
    image_attribute_data['image_time_data'] = np.zeros((n_images,))

    # Check which is last algo in list
    # This is indication of the data are detector images or spectra.
    # It also indicates if output is integer and binned
    last_algo = algo_list[-1]
    if is_detector_image(last_algo):

        is_binned = False
        if 'Binning' in algo_list:
            is_binned = True
        create_detector_image_output(output, dimensions, image_attribute_data, is_binned=is_binned)

    else:
        create_observation_data_output(output, dimensions, image_attribute_data)

    output_datasets.add_container('final',output)

# Inbetween output

    for algo_name in algo_list:

        # These are the inbetween output files
        # Maybe add some switch if we want them or not.
#        algo_output = input_data.get_dataset('im_algo_output', c_name='config', group='io')
        if kind == 'l1a':
            algo_output = input_data.get_dataset('im_algo_output', c_name='config', group='io')
        else:
            algo_output = input_data.get_dataset('l1b_algo_output', c_name='config', group='io')

        algo_file = algo_output.format(algo_name=algo_name)
        output_algo = dn.DataNetCDF(log, algo_file)
        add_main_attributes(output_algo, attributes)

        if is_detector_image(algo_name):
            is_binned = False
            if algo_name in ['Binning','ADC']:
                is_binned = True
            create_detector_image_output(output_algo, dimensions, image_attribute_data,
                                         is_binned=is_binned, in_between=True)

        else:
            create_observation_data_output(output_algo, dimensions, image_attribute_data,
                                           in_between=True)

        output_datasets.add_container(algo_name,output_algo)

    return output_datasets


def instrument_model(config, attributes):
    """
        Instrument model.
        Apply effects to create L0 data
    """

    start_time_im = time.perf_counter()

    # Get the input data into list of data containers
    input_data = get_input_data(config)

    # Get the list of algorithms from proctable file
    algo_list_input = input_data.get_dataset('algo_list', c_name='config', group='proctable')
    algo_list = input_data.get_dataset(algo_list_input, c_name='proctable')
    log.info(f"NOW algo_list: {algo_list}")

    # Get the input data. Note: possibly this needs to be updated when SGM is updated.
    radiance_data = input_data.get_dataset('radiance', c_name='measurement', kind='variable')

    #Note: this is the along track dimension
#    n_images = radiance_data.shape[0]
    # requested number of images:
    image_start = 0
    image_end = radiance_data.shape[0] - 1
    if 'image_start' in config:
        image_start = config['image_start']
    if 'image_end' in config:
        image_end = config['image_end']
    # image_end is up and including
    n_images = image_end+1 - image_start

    dim_alt = n_images
    dim_spat = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
    dim_spec = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
    dim_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')

    # Find binned row dimension
    bin_id = input_data.get_dataset('binning_table_id', c_name='config',
                                    group='detector', kind='variable')
    table = f"Table_{bin_id}"
    nr_binned_pixels = input_data.get_dataset('bins', c_name='binning',
                                              group=table, kind='dimension')
    binned_rows = int(nr_binned_pixels/dim_spec)

    log.info(f"Found dim_spec: {dim_spec} and dim_spat: {dim_spat} and dim_act: {dim_act}, and binned_rows: {binned_rows}")
    dimensions = {'dim_act':dim_act,'dim_spec': dim_spec, 'dim_spat': dim_spat,
                  'dim_alt': dim_alt, 'dim_binned_rows': binned_rows}

    # Get the output data into list of data containers
    output_datasets = initialize_output('l1a', input_data, algo_list, dimensions, attributes)

    # Loop over images
    for img in range(image_start, image_end+1):

        print(f"Processing image: {img}")
        start_time_image = time.perf_counter()

        image = radiance_data[img,:,:]
        # Have to temporary store the image somewhere for acces later on in algos.
        input_data.add_container('work', {'image': image})

        for algo_name in algo_list:
            start_time_algo = time.perf_counter()

            log.debug(f"IMG: {img} running algo {algo_name}")

            algo = globals()[algo_name](log)
            algo.check_input(input_data)

            log.debug(f"BEFORE RUNNING ALGO SHAPE IMAGE: {image.shape}")

            algo.execute( input_data)
            image = algo.get_data()

            log.debug(f"AFTER RUNNING ALGO SHAPE IMAGE: {image.shape}")

            log.debug(f"Updating measurement dataset for container {algo_name} with image nr {img}")
            output_datasets.update_dataset('measurement', data=image, c_name=algo_name, img=img)

            input_data.update_dataset('image', 'work',image)
            this_dataset = output_datasets.get_dataset('measurement', c_name=algo_name,
                                                       kind='variable')

            log.debug(f"{algo_name} SHAPE : {this_dataset.shape}")

            end_time_algo = time.perf_counter()
            log.info(f"Processing algorithm {algo_name} for image {img} took {(end_time_algo - start_time_algo):.6f}s")

        end_time_image = time.perf_counter()
        log.info(f"Processing image {img} took {(end_time_image - start_time_image):.6f}seconds")

        for algo_name in algo_list:
            algo_dataset = output_datasets.get_dataset('measurement', c_name=algo_name,
                                                       kind='variable')
            if is_detector_image(algo_name):
                output_datasets.update_dataset('detector_image', data=algo_dataset,
                                               c_name=algo_name, group='science_data')
            else:
                output_datasets.update_dataset('i', data=algo_dataset, c_name=algo_name,
                                               group='observation_data')

    final_dataset = output_datasets.get_dataset('measurement', c_name=algo_list[-1], kind='variable')
    if is_detector_image(algo_list[-1]):
        output_datasets.update_dataset('detector_image', data=final_dataset, c_name='final',
                                       group='science_data')
    else:
        output_datasets.update_dataset('i', data=final_dataset, c_name='final',
                                       group='observation_data')

    output_datasets.write()
    end_time_im = time.perf_counter()
    log.info(f"Instrument Model Processing took {(end_time_im - start_time_im):.6f}seconds")

    log.info("Succes")

if __name__ == '__main__' :

    # Get configuration info
    cfg_file = sys.argv[1]
    full_config = Utils.get_config(cfg_file)

    im_only_config = get_im_config(full_config)

    # Get information (like git hash and config file name and version (if available)
    # that will be added to the output file as attributes
    main_attributes = Utils.get_main_attributes(im_only_config, config_attributes_name='IM_configuration')

    # Also get some other attributes?
    im_config = get_specific_config(full_config, 'IM')
    attribute_dict = add_module_specific_attributes(im_config, main_attributes, 'im')

    # Get started
    instrument_model(im_config, attribute_dict)
