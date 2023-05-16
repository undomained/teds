import numpy as np
import sys

class surface_prop:
    """
    # The surface_prop class collects methods to
    # calculate surface properties
    #
    # CONTAINS
    # method __init__(self,wave)
    # method get_albedo_poly(self, albedo_coeffcients)
    # method get_albedo_CUSTOM(filename,sfname)
    # method get_albedo_ECOSTRESS(filename)
    """
    ###########################################################
    def __init__(self, wave):
        """
        # init class
        #
        # arguments: 
        #            wave: array of wavelengths [wavelength] [nm]
        """
        wave_max  = np.amax(wave)
        wave_min  = np.amin(wave)
        wave_mean = np.mean(wave)
        spec      = (wave-wave_mean)/(wave_max-wave_min)
        self.spec = spec
        self.alb  = np.zeros_like(wave)
        
    ###########################################################
    def get_albedo_poly(self, alb_coeff):
        """
         Parameters
        ----------
        spec : spectral array (can be any spectral quantity)
        alb  : list of albedo coefficients

        Returns
        -------
        albedo polynomial dependent on wavelenght

        """
        albedo = np.zeros(len(self.spec))  
        for i in range(0,len(alb_coeff)):
            albedo=albedo + alb_coeff[i]*(self.spec)**(i) 
            
        self.alb=albedo
        
    ###########################################################
    def get_albedo_CUSTOM(self,filename,sfname):
        """
        # Read albedo from custom database. This is generic typical 
        # data. For a comprehensive albedo database, have a look 
        # at: https://speclib.jpl.nasa.gov/
        #    
        # arguments: 
        #            filename: file with albedo database
        #            sfname: name of surface type 
        #                    [sand,soil,snow,vegetation,water]
        # returns:  
        #            alb: albedo array interpolated to wavelength [wavelength]
        """    
        #check whether input is in range
        sftypes=['sand','soil','snow','vegetation','water']
        while True:
            if os.path.exists(filename) and sfname in sftypes:
               break
            else:
               print("ERROR! surface_prop.get_albedo_CUSTOM: input out of range.")
               raise StopExecution
    
        # Read data from file
        raw = np.genfromtxt(filename,skip_header=15) # read file into numpy array, skip 15 header lines
        
        # Index of surface type
        isf = sftypes.index(sfname)
        
        # Interpolate albedo to wavelength array
        wave_org = raw[:,0]
        alb_org = raw[:,isf+1]
        self.alb = np.interp(self.wave,wave_org,alb_org)

    ###########################################################
    def get_albedo_ECOSTRESS(self,filename):
        """
        # Read albedo from ECOSTRESS database
        # at: https://speclib.jpl.nasa.gov/
        #
        # arguments: 
        #            filename: file with albedo database
        # returns:  
        #            alb: albedo array interpolated to wavelength [wavelength]
        """ 
        # Check whether input is in range
        while True:
            if os.path.exists(filename):
               break
            else:
               print("ERROR! surface_prop.get_albedo_ECOSTRESS: input out of range.")
               raise StopExecution
            
        # Read data from file
        raw = np.genfromtxt(filename,skip_header=21,unpack=True) 
    
        wv_in = np.array([a*1E3 for a in raw[0,:]]) # wavelength [nm]
        alb_in = np.array([a/1E2 for a in raw[1,:]]) # albedo [0...1]
        # Check if wavelength in ascending order. If not, flip arrays.
        if wv_in[0] > wv_in[-1]:
           # Interpolate albedo to wavelength array
           self.alb = np.interp(self.wave,np.flip(wv_in),np.flip(alb_in))
        else:
           self.alb = np.interp(self.wave,wv_in,alb_in)
