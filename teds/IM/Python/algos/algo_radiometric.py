import numpy as np
#from scipy.interpolate import CubicSpline

from teds.IM.Python.algos.algo_base import Algorithm

class Radiometric(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Radiometric algoritm specific code

    """

    def __init__(self, logger, algo_name="Radiometric"):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None
        self._stdev = None

    def check_input(self, input_data):
        """
            Check input data
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        # TODO: What would be a usefull check?

    def execute(self, input_data, kind='IM'):
        """
            Execute the algorithm
        """

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        self._stdev = stdev

        radiometric = input_data.get_dataset('radiometric', c_name='ckd', group='radiometric', kind='variable')

        exposure_time = input_data.get_dataset('exposure_time', c_name='config', group='detector', kind='variable')

        factor = (1.0/(exposure_time)) * radiometric
        if kind == 'IM':
            new_image = np.divide(image,factor)
        else:
            new_image = np.multiply(image,factor)
            new_stdev = np.multiply(stdev,factor)
            self._stdev = new_stdev

        self._data = new_image

        return




