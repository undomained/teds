import numpy as np


def generate(ncc):
    """
    Radiometric Calibration Constant
    """
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']

    dispersion = ncc.get_var("radiometric/dispersion") # in um/nm
    transmission = ncc.get_var("radiometric/transmission") # in fraction
    etendue =  cfg["etendue"]
    radiometric = 1 / (dispersion * transmission * etendue)  # nm um-1 mm-2 sr-2 
    
    #rcc = np.ones(ncc.get_shape(dims)) * cfg['rcc']
    attr = {
        'long_name': 'Radiometric calibration constant',
        'comment' : 'Transmission(wl) * Dispersion(wl) * Solid Angle',
        'comment2' : "wavelength independend single valued: FOV_act * FOV_ALT * pupil_act * pupil_alt * altitude",
        'units' : ""
    }
   
    ncc.create_var_auto(dims, radiometric, attr, 'f8')

