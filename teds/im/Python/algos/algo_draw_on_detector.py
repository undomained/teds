import numpy as np
from scipy.interpolate import CubicSpline

from teds.im.Python.algos.algo_base import Algorithm

class Draw_On_Detector(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Draw_On_Detector algoritm specific code

    """

    def __init__(self, logger, algo_name="Draw_On_Detector"):
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
            Execute the algorithm
            Cubic Spline is used to 'translate' across track to detector rows
        """

        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        n_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')
        n_row = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
        n_col = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
        row_indices = input_data.get_dataset('row_index', c_name='ckd', group='swath', kind='variable')
        # Note: is exected to be increasing
        if row_indices[-1,0] < row_indices[0,0]:
            # Need to reverse
            # Maybe also need to reverse image???????????????
            row_indices = np.flip(row_indices,axis=0)

        det_ispat = np.arange(n_row)

        detector_data = np.zeros((n_row, n_col))

        #When we CAN NOT assume row_indices to be the same for all ispec:
        for i_spec in range(n_col):
            x = row_indices[:,i_spec]
            y = image[:,i_spec]
            cs = CubicSpline(x,y, bc_type='natural')
            xs = det_ispat
            detector_data[:,i_spec] = cs(xs)

#        #When we can assume row_indices to be the same for all ispec:
#        cs = CubicSpline(row_indices[:,0], image,axis=0)
#        y = cs(det_ispat)
#        detector_data = y

        self._data = detector_data

        return




