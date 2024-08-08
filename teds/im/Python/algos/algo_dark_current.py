import copy

from teds.im.Python.algos.algo_base import Algorithm

class Dark_Current(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Dark_Current algoritm specific code

    """

    def __init__(self, logger, algo_name='Dark_Current'):
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
            kind IM: Add dark current to the data.
            kind L1B: Subtract dark current of the data.
            Only if enabled.
            Only for good pixels
        """
        self._logger.debug(f"Execute code from {self._algo_name} class")

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        if stdev is not None:
            self._stdev = stdev
        enabled = input_data.get_dataset('enabled', c_name='config', group='dark')
        if not enabled:
            # Algorithm will not be run
            self._logger.info(f"Algorithm {self._name} will not ne run because enabled is set to {enabled} in configuration file")
            return

        dark_current = input_data.get_dataset('current', c_name='ckd', group='dark', kind='variable')
        exposure_time = input_data.get_dataset('exposure_time', c_name='config', group='detector', kind='variable')

        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')
        new_image = copy.deepcopy(image)
        good_pixels = pixel_mask==0
        if kind == 'IM':
            new_image[good_pixels] += exposure_time * dark_current[good_pixels]
        else:
            new_image[good_pixels] -= exposure_time * dark_current[good_pixels]
        self._data = new_image

        return 

