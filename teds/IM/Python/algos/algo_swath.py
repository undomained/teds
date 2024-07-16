import numpy as np
#from scipy.interpolate import CubicSpline

from teds.IM.Python.algos.algo_base import Algorithm

class Swath(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with Swath algoritm specific code

    """

    def __init__(self, logger, algo_name="Swath"):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None
        self._stdev = None
        self._mask = None

    def check_input(self, input_data):
        """
            Check input data
        """
        self._logger.debug(f"Check INPUT from {self._algo_name} class")
        # TODO: What would be a usefull check?

    def calculate_pixel_indices_and_weights(self,input_data, spectrum_width = 1):
        """
           Translated directly from C++ code. 
           Calculating pixel indices and weights
        """
        spectrum_width_half = spectrum_width/2
        n_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')
        n_col = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
        row_indices = input_data.get_dataset('row_index', c_name='ckd', group='swath', kind='variable')

        n_indices = spectrum_width + 1
        swath_pix_indices = np.zeros((n_act, n_col*n_indices))
        swath_weights = np.zeros((n_act, n_col*n_indices))

        # Based on the argument spectrum_width, each spectral element
        # spans N rows on the detector. The weight of the first and last
        # rows are determined by where exactly the edge of the spectrum
        # falls. All middle rows have the same constant weight (1 before
        # normalization).


        for i_act in range(n_act):
            indices = swath_pix_indices[i_act,:]
            weights = swath_weights[i_act,:]
            for i in range(n_col):
                row_first_d = row_indices[i_act][i] + 0.5 - spectrum_width_half
                row_last_d = row_indices[i_act][i] + 0.5 + spectrum_width_half
                row_first_i =  int(row_first_d)
                row_last_i  = int(row_last_d)

                indices[i * n_indices] = row_first_i * n_col + i
                indices[(i + 1) * n_indices - 1] = row_last_i * n_col + i
                weights[i * n_indices] = np.ceil(row_first_d) - row_first_d
                weights[(i + 1) * n_indices - 1] = row_last_d - np.floor(row_last_d)
                # Get weights for inner indices
                for i_ind in range( 1, n_indices):
                    idx = i * n_indices + i_ind
                    indices[idx] = (row_first_i + i_ind) * n_col + i;
                    weights[idx] = 1.0
                # weights outside
                idx = i * n_indices
                if (indices[idx] == indices[idx + 1]):
                    weights[idx + 1] = 0.0
                if (indices[idx + n_indices - 1] == indices[idx + n_indices - 2] and n_indices > 2):
                    weights[idx + n_indices - 2] = 0.0;

                weights_norm = 0
                for i_ind in range(n_indices):
                    weights_norm += weights[i * n_indices + i_ind]

                for i_ind in range(n_indices):
                    weights[i * n_indices + i_ind] /= weights_norm;

            swath_pix_indices[i_act,:] = indices
            swath_weights[i_act,:] = weights

        return swath_pix_indices, swath_weights


    def execute(self, input_data, kind='L1AL1B', spectrum_width=1):
        """
            Execute the algorithm
        """

        swath_pix_indices, swath_weights = self.calculate_pixel_indices_and_weights(input_data)

        image = input_data.get_dataset('image', c_name='work')
        stdev = input_data.get_dataset('stdev', c_name='work')
        self._data = image
        self._stdev = stdev

        n_act = input_data.get_dataset('across_track', c_name='ckd', kind='dimension')
        n_row = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
        n_col = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
        row_indices = input_data.get_dataset('row_index', c_name='ckd', group='swath', kind='variable')
        pixel_mask = input_data.get_dataset('pixel_mask', c_name='ckd', kind='variable')

        spectrum_width_half = spectrum_width/2
        n_indices = spectrum_width + 1

        new_image = np.zeros((n_act, n_col))
        new_stdev = np.zeros((n_act, n_col))
        new_mask = np.zeros((n_act, n_col))

        # Like this? Or use flatten?
        image_reshaped = np.reshape(image, (n_row*n_col,))
        stdev_reshaped = np.reshape(stdev, (n_row*n_col,))
        pixel_mask_reshaped = np.reshape(pixel_mask, (n_row*n_col,))

        for i_act in range(n_act):
            for i in range(n_col):
                total_weight = 0
                for i_ind in range(n_indices):
                    idx = i * n_indices + i_ind
                    ipix =  int(swath_pix_indices[i_act,idx])
                    # Note ipix is index in 2D array. n_row*n_col
                    # reshaping is needed
    
                    # Protection against ipix values that are outside the pixel_mask dimension
                    if ((ipix > 0) and (ipix <  pixel_mask_reshaped.shape[0])):
                        if not pixel_mask_reshaped[ipix]:
                            w  = swath_weights[i_act,idx]
                            new_image[i_act, i] += image_reshaped[ipix] * w;
                            n = stdev_reshaped[ipix] * w 
                            new_stdev[i_act,i] += n * n
                            total_weight += w
                if (total_weight > 0.0):
                    new_image[i_act,i] /= total_weight
                    new_stdev[i_act,i] = np.sqrt(new_stdev[i_act,i] / (total_weight * total_weight))
                else:
                    new_mask[i_act,i] = True

        self._data = new_image
        self._stdev = new_stdev
        self._mask = new_mask

