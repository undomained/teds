import copy

from teds.im.Python.algos.algo_base import Algorithm
from teds import log

class Dark_Offset(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Offset algoritm specific code

    """

    def __init__(self, algo_name='Dark_Offset'):
        
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
            kind IM: Add dark offset to the data.
            kind L1B: Subtract dark offset from the data.
            Only if enabled.
            Only for good pixels
        """
        log.debug(f"Execute code from {self._algo_name} class")

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        if stdev is not None:
            self._stdev = stdev
        enabled = input_data.get_dataset('enabled', c_name='config', group='dark')
        if not enabled:
            # Algorithm will not be run
            log.info(f"Algorithm {self._name} will not ne run because enabled is set to {enabled} in configuration file")
            return
        dark_offset = input_data.get_dataset('offset', c_name='ckd', group='dark', kind='variable')
        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')
        new_image = copy.deepcopy(image)
        good_pixels = pixel_mask==0
        if kind == 'IM':
            new_image[good_pixels] += dark_offset[good_pixels]
        else:
            new_image[good_pixels] -= dark_offset[good_pixels]
        self._data = new_image

        return 

