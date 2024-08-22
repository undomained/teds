import numpy as np

def generate(ncc):
    """
    Radiometric Calibration Constant
    """
    cfg = ncc.cfg
    dims = ['detector_row', 'detector_column']
    rcc = np.ones(ncc.get_shape(dims)) * cfg['rcc']
    attr = {
        'long_name': 'Radiometric calibration constant',
        'comment' : 'Single value, will be replaced with \
            wavelength dependent dispersion and transmission'
    }
    
    ncc.create_var_auto(dims, rcc, attr, 'f8')
