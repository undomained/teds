import sys
import numpy as np
from teds.im.Python.algos.algo_base import Algorithm

class PRNU(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with PRNU algoritm specific code
    """

    def __init__(self, logger, algo_name="PRNU"):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None
        self._stdev = None


    def check_input(self, input_data):
        """
            Check on input data.
            Dimensions of data and PRNU should be the same
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        # TODO do we want to call it image or data???
        # data is more general name
        image = input_data.get_dataset('image', c_name='work')
        if image is None:
            error_message = f"No input data given for algorithm {self._algo_name}. Exiting."
            self._logger.error(error_message)
            sys.exit(error_message)
        prnu = input_data.get_dataset('prnu', c_name='ckd', group='prnu', kind='variable')
        if prnu is None:
            error_message = f"No prnu ckd data given for algorithm {self._algo_name}. Exiting."
            self._logger.error(error_message)
            sys.exit(error_message)

        data = image
        prnu_ckd = prnu

        self._logger.debug(f"DATA shape: {data.shape} AND CKD shape: {prnu_ckd.shape}")
        # SHAPES SHOULD BE THE SAME
        if data.shape != prnu_ckd.shape:
            error_message = f"PRNU shape {prnu_ckd.shape} and data shape {data.shape} not the same! Exiting!"
            self._logger.error(error_message)
            sys.exit(error_message)
        return

    def execute(self, input_data, kind='IM'):
        """
            Perform the correction.
            Only if enabled.
            Only for good pixels
        """

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        if stdev is not None:
            self._stdev = stdev
        enabled = input_data.get_dataset('enabled', c_name='config', group='prnu')
        if not enabled:
            # Algorithm will not be run
            self._logger.info(f"Algorithm {self._name} will not ne run because enabled is set to {enabled} in configuration file")
            return
        prnu_ckd = input_data.get_dataset('prnu', c_name='ckd', group='prnu', kind='variable')
        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')
        self._logger.debug(f"Execute code from {self._algo_name} class")
        if kind == 'IM':
            new_image = np.multiply(image,prnu_ckd, where=pixel_mask==0)
        else:
            new_image = np.divide(image,prnu_ckd, where=pixel_mask==0)
            new_stdev = np.divide(stdev,prnu_ckd, where=pixel_mask==0)
            self._stdev = new_stdev
        self._data = new_image

        return

