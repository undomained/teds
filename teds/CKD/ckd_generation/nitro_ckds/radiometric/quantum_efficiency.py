"""
Quantum Efficiency
Wavelength dependend -> Uses Wavelength Map
"""

from scipy.interpolate import interp1d

def generate(ncc):
    dims = ['detector_row', 'detector_column']

    # data, source: TANGO Nitro SNR budget.xlsx (TNO)
    wl = [405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490]
    qe = [0.54, 0.56, 0.58, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.685, 0.70, 0.7, 0.705, 0.71, 0.71]
    f = interp1d(wl,qe, fill_value = 'extrapolate')

    # get or calculate wavelength map
    wl_map = ncc.get_var('spectral/wave_map')
    QE = f(wl_map)
    attr = {
        'comment' : 'Generated using wavelength map and QE data from TNO'
    }
    
    ncc.create_var_auto(dims, QE, attr, 'f8')






