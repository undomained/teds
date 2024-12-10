import numpy as np


def generate(ncc):
    """
    Radiometric Calibration Constant
    """
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']

    dispersion = ncc.get_var("radiometric/dispersion") # in um/nm
    transmission = ncc.get_var("radiometric/transmission") # in fraction
    etendue =  cfg["etendue"] * 10**(-6) # m2 sr 
    dl = 1/dispersion * cfg["detector_pixel_size"]*10**3 # dlambda in nm/px
    gain = cfg["gain"]

    rad_factor = dl * transmission * etendue / gain # nm/px m2 sr
    radiometric = 1 / rad_factor

    print("rcc = {:.2e}, new radiometric = {:.2e}".format(cfg['rcc'], np.mean(radiometric)))

    # TODO: WORK IN PROGRESS

    # old radiometric constant
    rcc = np.ones(ncc.get_shape(dims)) * cfg['rcc']


    attr = {
        'long_name': 'Radiometric calibration constant',
        'comment' : 'Transmission(wl) * Dispersion(wl) * Solid Angle / Gain',
        'comment2' : "Still constant is used, will be wavelength dependend",
        'units' : ""
    }
   
   # for now use the old constant
    ncc.create_var_auto(dims, rcc, attr, 'f8')

