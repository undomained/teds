import numpy as np
from teds import log
from teds.im.Python.algos.algo_base import Algorithm

class ISRF(Algorithm):
    """
        Sub class of base class Algorithm
        Base class methods check_input and execute are overwritten with ISRF algoritm specific code
        Extra Methods:
        - apply_lin_interpol(self, input_data)
        - apply_isrf(self,input_data)

    """

    def __init__(self, algo_name="ISRF"):
        
        self._algo_name = algo_name
        self._data = None
        self._do_wldep_isrf = False
        log.debug("INIT ISRF module")


    def check_input(self, input_data):
        """
            Check input data
        """
        log.debug(f"Check INPUT from {self._algo_name} class")

        # Check if input image and ISRF CKD have the same wavelength grids
        image = input_data.get_dataset('image', c_name = 'work')
        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
        n_act = image.shape[0]
        wavelength = np.tile(wavelength, (n_act,1)) # tile wavelengths to act_pos

        # 3D isrf data, one ISRF per wavelength. 
        # Wavelength grid should be equal to wavelength grid of input radiance.
        isrf_args = {"c_name" : 'ckd', "group" : 'radiometric', "kind" : 'variable'}
        isrf    = input_data.get_dataset('isrf', **isrf_args)
        isrf_wl = input_data.get_dataset('isrf_wavelengths', **isrf_args)
        
        # Set correct 
        dimsizes = np.array(isrf.shape)
        dim0 = np.flatnonzero(dimsizes == wavelength.shape[0])
        dim1 = np.flatnonzero(dimsizes == wavelength.shape[1])
        
        if not (len(dim0) or not len(dim1)): # no matches in dimensions size
            msg = f"ISRF CKD wavelength grid does not match input wavelengths, using constant ISRF"
            log.warning(msg)
            self._do_wldep_isrf = False
        else: # put dimensions in correct order
            isrf = np.transpose(isrf, (dim0[0], dim1[0], 2))
            isrf_wl = np.transpose(isrf_wl, (dim0[0], dim1[0], 2))
            input_data.update_dataset('isrf', 'ckd', isrf, group = 'radiometric')
            input_data.update_dataset('isrf_wavelengths', 'ckd', isrf_wl, group = 'radiometric')

            center_ix = int(isrf_wl.shape[2]/2)
            isrf_wl_center = isrf_wl[:,:, center_ix] # get center wavelength of ISRF shapes

            # shape wavelength grids are same, now check the values
            diff_wl = 0
            for i_act in range(wavelength.shape[1]):
                # find difference in wavelength up to 3 decimals
                diff_per_act = np.round(np.abs(isrf_wl_center[:,i_act] - wavelength[:,i_act]), 3)
                # add sum over the entire act pos
                diff_wl += np.sum(diff_per_act)
            if not diff_wl:
                self._do_wldep_isrf = True
            else:
                msg = f"ISRF CKD wavelength values do not match input wavelengths, using constant ISRF"
                log.warning(msg)
                self._do_wldep_isrf = False
    
    def apply_lin_interpol(self, input_data):
        """
            Linearly interpolate the line-by-line spectra onto 
            the CKD wavelength grids.
        """
        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        # Wavelength that belongs to the image.
        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')

        # Wavelength map from ckd
        wavemap = input_data.get_dataset('wavelength', c_name='ckd', group='spectral', kind='variable')

        interp_data = np.zeros(wavemap.shape)
        for i_act in range(image.shape[0]):
            interp_data[i_act,:] = np.interp(wavemap[i_act,:], wavelength, image[i_act,:])

        return interp_data
                         

    def apply_isrf_wldep(self, input_data):
        """
            Convolve the radiances with an ISRF from the line-by-line grid
            onto the target grids taken from the CKD.
            Using wavelength dependent ISRF shape from CKD
        """
        
        log.debug('Applying wavelength dependent ISRF')

        image = input_data.get_dataset('image', c_name = 'work')
        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
        wavelength = np.tile(wavelength, (image.shape[0],1)) # tile wavelengths to act_pos

        # Wavelength grid should be equal to wavelength grid of input radiance.
        isrf_args = {"c_name" : 'ckd', "group" : 'radiometric', "kind" : 'variable'}
        isrf = input_data.get_dataset('isrf', **isrf_args)
        isrf_wl = input_data.get_dataset('isrf_wavelengths', **isrf_args)

        n_actpos = isrf_wl.shape[0]
        n_wavelengths = isrf_wl.shape[1]
        n_samples_isrf = isrf_wl.shape[2]
        center_ix = int(n_samples_isrf/2)

        image_conv = np.zeros(image.shape, dtype = float)
        for i_act in range(n_actpos):
            for i_wl in range(n_wavelengths):
                # left/right part of first/last isrf samples is outside wl range
                i_isrf_0 = np.max([center_ix - i_wl, 0]) 
                i_isrf_1 = np.min([n_wavelengths + center_ix - i_wl - 1, n_samples_isrf - 1])
                i_img_0 = np.max([i_wl - center_ix, 0])
                i_img_1 = np.min([i_wl + center_ix, n_wavelengths - 1])

                this_isrf = isrf[i_act, i_wl, i_isrf_0:i_isrf_1]
                weights =  this_isrf / np.sum(this_isrf) 

                # apply convolution
                image_conv[i_act, i_wl] = np.sum(image[i_act, i_img_0:i_img_1] * weights)

        return image_conv



    def apply_isrf_fancy(self,input_data):
        """
            Convolve the radiances with an ISRF from the line-by-line grid
            onto the target grids taken from the CKD.
            Using 1 FWHM given by configuration file
            Apply Gauss function on each element of array.
        """

        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        # Wavelength that belongs to the image.
        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
        log.debug(f"WAVELENGTH SHAPE At THE MOMENT: {wavelength.shape}")
        n_act = image.shape[0]
        wavelength = np.tile(wavelength, (n_act,1))
        log.debug(f"NEW SHAPE: {wavelength.shape}")

        # Wavelength map from ckd
        wavemap = input_data.get_dataset('wavelength', c_name='ckd', group='spectral', kind='variable')
        log.debug(f"WAVEMAP SHAPE: {wavemap.shape}")

        fwhm_gauss = input_data.get_dataset('fwhm_gauss', c_name='config', group='isrf')
        log.debug(f"input fwhm: {fwhm_gauss}")

        # number of detector columns obtained from ckd
        n_det_cols = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')

        # Create image with new dimensions
        image_conv = np.zeros((image.shape[0], n_det_cols))

        # calculate factors to use in gauss calculation
        sigma = fwhm_gauss/(2.*np.sqrt(2.*np.log(2.)))
        sigma_inv = 1.0/(2.*sigma*sigma)
        norm_inv = 1.0/(sigma * np.sqrt(2*np.pi))

        # Create lookup table with weights
        fwhm = np.ones(wavelength.shape)*fwhm_gauss
        # a la C++ code use a limit for which the gauss weight is set to zero.
        wave_limit = 3.0 * fwhm
        def apply_gauss(x):
            return (norm_inv*np.exp(-(x-mu)*(x-mu)*sigma_inv))

        delta_lambdas = np.zeros(wavelength.shape)
        for i_act in range(wavelength.shape[0]):
            values = wavelength[i_act,:]
            delta_lambda = 0.5*(values[2:] - values[:-2])
            # add first and last value to create correct shape
            delta_lambda = np.insert(delta_lambda, 0, (values[1]-values[0]))
            delta_lambda = np.append(delta_lambda, values[-1] - values[-2])
            delta_lambdas[i_act,:] = delta_lambda

        for i_det in range(n_det_cols):
            mu = np.column_stack([wavemap[:,i_det].tolist()]*wavelength.shape[1])
            #Apply Gauss function on each element of array.
            gauss_weights = apply_gauss(wavelength)
            too_far = abs(wavelength - mu) > wave_limit
            gauss_weights[too_far] = 0.
            weights = np.multiply(gauss_weights, delta_lambdas)
            conv_signal = np.multiply(weights, image)
            image_conv[:,i_det] = np.sum(conv_signal, axis=1)

        return image_conv

    def apply_isrf(self,input_data):
        """
            Convolve the radiances with an ISRF from the line-by-line grid
            onto the target grids taken from the CKD.
            Using 1 FWHM given by configuration file
            Try to be a bit faster
        """

        image = input_data.get_dataset('image', c_name='work')
        self._data = image

        # Wavelength that belongs to the image.
        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
        log.debug(f"WAVELENGTH SHAPE At THE MOMENT: {wavelength.shape}")
        n_act = image.shape[0]
        wavelength = np.tile(wavelength, (n_act,1))
        log.debug(f"NEW SHAPE: {wavelength.shape}")

        # Wavelength map from ckd
        wavemap = input_data.get_dataset('wavelength', c_name='ckd', group='spectral', kind='variable')
        #log.debug(f"WAVEMAP SHAPE: {wavemap.shape}")

        fwhm_gauss = input_data.get_dataset('fwhm_gauss', c_name='config', group='isrf')
        log.debug(f"input fwhm: {fwhm_gauss}")

        # number of detector columns obtained from ckd
        n_det_cols = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')

        # Create image with new dimensions
        image_conv = np.zeros((image.shape))

        # calculate factors to use in gauss calculation
        sigma = fwhm_gauss/(2.*np.sqrt(2.*np.log(2.)))
        sigma_inv = 1.0/(2.*sigma*sigma)
        norm_inv = 1.0/(sigma * np.sqrt(2*np.pi))

        # Create lookup table with weights
        fwhm = np.ones(wavelength.shape)*fwhm_gauss
        wave_limit = 3.0 * fwhm

        delta_lambdas = np.zeros(wavelength.shape)
        for i_act in range(wavelength.shape[0]):
            values = wavelength[i_act,:]
            delta_lambda = 0.5*(values[2:] - values[:-2])
            # add first and last value to create correct shape
            delta_lambda = np.insert(delta_lambda, 0, (values[1]-values[0]))
            delta_lambda = np.append(delta_lambda, values[-1] - values[-2])
            delta_lambdas[i_act,:] = delta_lambda

        for i_det in range(n_det_cols):
            mu = np.column_stack([wavemap[:,i_det].tolist()]*wavelength.shape[1])
            wave_diff = wavelength - mu
            too_large = abs(wave_diff) > wave_limit

            wave_diff[abs(wave_diff)>wave_limit] = 0.

            # Using where to speed up a little bit
            x2 = np.square(wave_diff, where=abs(wave_diff)>0.)
            gauss_weights = norm_inv*np.exp(-sigma_inv*x2, where=x2>0.)

            gauss_weights[too_large] = 0.
            weights  = np.multiply(gauss_weights, delta_lambdas)
            conv_signal = np.multiply(weights, image)
            image_conv[:,i_det] = np.sum(conv_signal, axis=1)

        return image_conv


    def execute(self, input_data):
        """
            Execute the algorithm.
            When enabled is set to 0 (False) only linear interpolation is performed.
            When enabled is set to 1 (True) wavelength concolution is performed
        """

        log.debug(f"Execute code from {self._algo_name} class")

        do_isrf = input_data.get_dataset('enabled', c_name='config', group='isrf')
        if do_isrf :
            if self._do_wldep_isrf:
                # apply wavelength dependent isrf if isrf-ckd and input matches
                self._data = self.apply_isrf_wldep(input_data) 
            else:
                # apply same isrf shape to every wavelength
                self._data = self.apply_isrf(input_data)
        else:
            self._data = self.apply_lin_interpol(input_data)
        
        return


