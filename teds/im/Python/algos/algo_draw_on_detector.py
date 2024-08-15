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

        wl_img = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
        
        n_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')
        n_row = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
        n_col = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
        
        
        wl_per_col = input_data.get_dataset('wavelength', c_name = 'ckd', group = 'spectral', kind = 'variable')
        row_indices = input_data.get_dataset('row_index', c_name='ckd', group='swath', kind='variable')
        
        # First convert wavelengths to columns
        lbl_in_cols = np.zeros((n_act, n_col))
        for i_act in range(n_act):
            spl_wl_to_col = CubicSpline(wl_per_col[i_act], np.arange(n_col))
            col_ix = spl_wl_to_col(wl_img) # decimal col indices of scene samples
            
            lbl = image[i_act] # line-by-line spectrum
            # Interpolate lbl to integer col indices
            spl_lbl_vs_col = CubicSpline(col_ix, lbl) 
            lbl_in_cols[i_act] = spl_lbl_vs_col(np.arange(n_col))
        
        # Now per column interpolate the rows
        det_img = np.zeros((n_row, n_col))
        for i_col in range(n_col):            
            lbl = lbl_in_cols[:, i_col] # scene per column
            spl_lbl_vs_row = CubicSpline(row_indices[:,i_col], lbl)
            # interpolate to integer rows
            det_img[:, i_col] = spl_lbl_vs_row(np.arange(n_row))
            
        self._data = det_img
        
        return



