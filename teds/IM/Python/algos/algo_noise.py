import copy
import numpy as np

from teds.IM.Python.algos.algo_base import Algorithm

class Noise(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Noise algoritm specific code

    """

    def __init__(self, logger, algo_name='Noise'):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None


    def check_input(self, input_data):
        """
            Check input data
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        # TODO: What would be a usefull check?

    def execute(self, input_data):
        """
            Add noise to the data.
            Only if enabled.
            Only for good pixels
        """
        self._logger.debug(f"Execute code from {self._algo_name} class")

        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        enabled = input_data.get_dataset('enabled', c_name='config', group='noise')

        if not enabled:
            # Algorithm will not be run
            self._logger.info(f"Algorithm {self._name} will not ne run because enabled is set to {enabled} in configuration file")
            return

        noise_n2 = input_data.get_dataset('n', c_name='ckd', group='noise', kind='variable')
        noise_g = input_data.get_dataset('g', c_name='ckd', group='noise', kind='variable')
        noise_seed = input_data.get_dataset('seed', c_name='config', group='noise', kind='variable')

        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')

        new_image = copy.deepcopy(image)

#        noise_value = np.sqrt(noise_n2 + image * noise_g)
        noise_value = np.sqrt(np.add(noise_n2, np.multiply(image,noise_g)))
        np.random.seed(noise_seed)
        random_normal_noise =  np.random.normal(0.0, noise_value)

        good_pixels = pixel_mask==0
        new_image[good_pixels] += random_normal_noise[good_pixels]
        self._data = new_image

        return 

