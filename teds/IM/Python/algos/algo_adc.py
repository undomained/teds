import sys
import numpy as np
from teds.IM.Python.algos.algo_base import Algorithm

class ADC(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with ADC algoritm specific code
    """

    def __init__(self, logger, algo_name="ADC"):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None


    def check_input(self, input_data):
        """
            Check on input data.
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        return

    def execute(self, input_data):
        """
            Perform the correction.
        """

        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        self._logger.debug(f"Execute code from {self._algo_name} class")
        new_image = np.rint(image).astype(int)
        self._data = new_image

        return

