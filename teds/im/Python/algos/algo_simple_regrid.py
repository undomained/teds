import numpy as np
from scipy.interpolate import CubicSpline

from teds.im.Python.algos.algo_base import Algorithm
from teds import log

class Simple_Regrid(Algorithm):

    def __init__(self, algo_name="Simple_Regrid"):
        
        self._algo_name = algo_name
        self._data = None
#        print("INIT Simple_Regrid module")

    def spatial_regrid(self, data, fov_ispat, det_ispat):

        # NOTE: THIS TAKES LONG TIME

        detector_data = np.zeros((det_ispat.shape[0], data.shape[1]))
#        for ispec in range(data.shape[1]):
#            x = fov_ispat[:,ispec]
#            y = data[:,ispec]
#            cs = CubicSpline(x,y)
#            xs = det_ispat
#            detector_data[:,ispec] = cs(xs)

        #fov_ispat is assumed to be the same for all ispec
        cs = CubicSpline(fov_ispat[:,0], data,axis=0)
        y = cs(det_ispat)
#        print(f"Y SHAPE HIER: {y.shape}")
        detector_data = y

        return detector_data

    def wave_regrid(self, data, source_grid, target_grid):
        interp_data = np.zeros(target_grid.shape)
        for ifov in range(data.shape[0]):
            interp_data[ifov,:] = np.interp(target_grid[ifov,:],source_grid,data[ifov,:])

        return interp_data


#    def check_input(self, image, ckd, wavelength):
    def check_input(self, input_data):
        print(f"Check INPUT from {self._algo_name} class")

#    def execute(self, image, ckd, wavelength):
    def execute(self, input_data):

        image = input_data.get_dataset('image', c_name='work')
        self._data = image
        wave_target = input_data.get_dataset('wave_target', c_name='ckd', group='WAVELENGTH', kind='variable')

        wavelength = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')
#        print(f"Execute code from {self._algo_name} class")
        # Preform simple regridding

        # wavelength regridding
        # For now get it from ckd. 
        # In Future get it from input_data
#        wave_target = ckd.get('wave_target', group='WAVELENGTH')
#        print(f"SHAPE WAVE_TARGET: {wave_target.shape}")
#        print(f"SHAPE SOURCE WAVELENGTH: {wavelength.shape}")

#        print(f"VALUES WAVE_TARGET: {wave_target}")
#        print(f"VALUES SOURCE WAVELENGTH: {wavelength}")

#        print(f"SHAPE IMAGE AT START: {image.shape}")

        image = self.wave_regrid(image, wavelength, wave_target)

#        print(f"SHAPE IMAGE AT END: {image.shape}")
        self._data = image

        # spatial regridding

        fov_ispat = input_data.get_dataset('fov_ispat', c_name='ckd', group='FIELD_OF_VIEW', kind='variable')
#        print(f"FOV_ISPAT SHAPE: {fov_ispat.shape}")
        # Note: is exected to be increasing
        if fov_ispat[-1,0] < fov_ispat[0,0]:
            # Need to reverse
            fov_ispat = np.flip(fov_ispat,axis=0)

#        print(f"FOV_ISPAT DATA[:,0]: {fov_ispat[:,0]}")
#        print(f"FOV_ISPAT DATA[:,100]: {fov_ispat[:,100]}")
#        print(f"FOV_ISPAT DATA[:,1500]: {fov_ispat[:,1500]}")

#        n_spat_det = ckd.get('spatial_detector_pixels', kind='dimension')
        n_spat_det = input_data.get_dataset('spatial_detector_pixels', c_name='ckd', kind='dimension')
#        print(f"N_SPAT_DET: {n_spat_det}")
        det_ispat = np.arange(n_spat_det)
#        print(f"DET_ISPAT SHAPE: {det_ispat.shape}")

#        print(f"IMAGE VALUES HIER: {image[:100]}")
        image = self.spatial_regrid( image, fov_ispat, det_ispat)
#        print(f"SHAPE IMAGE AT END2: {image.shape}")
#        print(f"IMAGE VALUES NOW: {image[:100]}")
        self._data = image

        return




