import copy
import numpy as np

from teds.im.Python.algos.algo_base import Algorithm
from teds import log

class Noise(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Noise algoritm specific code

    """

    def __init__(self, algo_name='Noise'):
        
        self._algo_name = algo_name
        self._data = None
        self._stdev = None


    def check_input(self, input_data):
        """
            Check input data
        """
        log.debug(f"Check INPUT from {self._algo_name} class")
        # TODO: What would be a usefull check?

    def execute(self, input_data, kind='IM'):
        """
            Add noise to the data.
            Only if enabled.
            Only for good pixels
        """
        log.debug(f"Execute code from {self._algo_name} class")

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        if stdev is not None:
            self._stdev = stdev

        enabled = input_data.get_dataset('enabled', c_name='config', group='noise')

        if not enabled:
            if kind != 'IM':
                stdev = np.ones(image.shape)
                self._stdev = stdev
            # Algorithm will not be run
            log.info(f"Algorithm {self._name} will not ne run because enabled is set to {enabled} in configuration file")
            return

        noise_n2 = input_data.get_dataset('n', c_name='ckd', group='noise', kind='variable')
        noise_g = input_data.get_dataset('g', c_name='ckd', group='noise', kind='variable')
        noise_seed = input_data.get_dataset('seed', c_name='config', group='noise', kind='variable')

        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')


        if kind == 'IM':
            new_image = copy.deepcopy(image)

#            noise_value = np.sqrt(noise_n2 + image * noise_g)
            noise_value = np.sqrt(np.add(noise_n2, np.multiply(image,noise_g)))
            np.random.seed(noise_seed)
            random_normal_noise =  np.random.normal(0.0, noise_value)

            good_pixels = pixel_mask==0
            new_image[good_pixels] += random_normal_noise[good_pixels]
            self._data = new_image
        else:
            # L1B determine noise level
            variance = np.add(noise_n2, np.multiply(image,noise_g))
            bad_pixels = variance < 0.
            pixel_mask[bad_pixels] = 1
            good_pixels = variance >= 0.

            # For L1B this info comes from l1a netcdf file
            nr_coadditions = input_data.get_dataset('nr_coadditions', c_name='measurement', group='image_attributes', kind='variable')[0]

            # For L1B this info comes from l1a netcdf file
            bin_id = input_data.get_dataset('binning_table', c_name='measurement', group='image_attributes', kind='variable')[0]

            table = f"Table_{int(bin_id)}" 
            count_table = input_data.get_dataset('count_table', c_name='binning', group=table, kind='variable')
            count_table = np.reshape(count_table, variance.shape)
            new_stdev = np.sqrt(np.divide(variance,np.multiply(count_table, nr_coadditions)))
            self._stdev = new_stdev
            
        return 

