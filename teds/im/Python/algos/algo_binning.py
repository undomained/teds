import sys
import numpy as np
from teds.im.Python.algos.algo_base import Algorithm

class Binning(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Binning algoritm specific code
    """

    def __init__(self, logger, algo_name="Binning"):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None


    def check_input(self, input_data):
        """
            Check on input data.
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        return

    def execute(self, input_data, other_data=None):
        """
            Perform the correction.
        """

        if other_data is not None:
            self._data = other_data
        else:
            image = input_data.get_dataset('image', c_name='work')
            self._data = image

        self._logger.debug(f"Execute code from {self._algo_name} class")

        bin_id = input_data.get_dataset('binning_table_id', c_name='config', group='detector')
        table = f"Table_{bin_id}" 
        print(f"Fetching binning table: {table}")
        count_table = input_data.get_dataset('count_table', c_name='binning', group=table, kind='variable')
        binning_table = input_data.get_dataset('binning_table', c_name='binning', group=table, kind='variable')
        binned_pixels = input_data.get_dataset('bins', c_name='binning', group=table, kind='dimension')
        det_cols = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
        det_rows = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
        binned_rows = int(binned_pixels/det_cols)

        binned_image = np.zeros((binned_rows, self._data.shape[1]))
        print(f"SHAPE OF BINNED IMAGE: {binned_image.shape}")

        # Note: nrs in binning-table are pixel numbers!
        for det_row in range(binning_table.shape[0]):
            pix_numbers = binning_table[det_row]
            binned_row = pix_numbers[0]
            # Taking a shortcut here. Assuming no binning in col direction.
            # Binned row number for all columns should be the same.
            if binned_row != 0:
                binned_row /= det_cols
            binned_image[int(binned_row),:] += self._data[det_row,:]

        self._data = binned_image

        return

