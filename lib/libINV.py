# libary of inversion routines

import numpy as np
import matplotlib.pyplot as plt
import sys

from end_to_end.lib import libRT
from end_to_end.lib import libSURF

###########################################################

def Gauss_Newton_iteration(retrieval_init, atm, optics, measurement, isrf, max_iter, chi2_lim):
    """
    Non-linear least square fit using Gauss-Newton iteration
    """

    k_dic = {}   # dictionary of the state vector, later inserted in nonscattered forward model
    x_dic = {}
    g_dic = {}
    dev = []  # correspsonding list for derivatives

    # initialization
    for key in retrieval_init['trace gases'].keys():
        if (key == 'CO2'):
            k_dic[key] = 'molec_07'
            dev.append('molec_07')
        if (key == 'CH4'):
            k_dic[key] = 'molec_32'
            dev.append('molec_32')
        if (key == 'H2O'):
            k_dic[key] = 'molec_01'
            dev.append('molec_01')
        scaling = retrieval_init['trace gases'][key]['scaling']
        Xcol = sum(atm.__getattribute__(key))/sum(atm.air)/scaling
        x_dic[key] = retrieval_init['trace gases'][key]['init']/Xcol   # prior scaling parameter
    m = 0
    while 'alb%d' % (m) in retrieval_init['surface'].keys():
        k_dic['alb%d' % (m)] = "alb%d" % (m)
        x_dic['alb%d' % (m)] = retrieval_init['surface']['alb%d' % (m)]
        m = m+1
    dev.append('alb')

    # start with iteration
    iteration = 0
    convergence = False
    Chi_sqrt = []

    surface = libSURF.surface_prop(retrieval_init['wavelength lbl'])

    convergence = False
    for iteration in range(retrieval_init['maximum iteration']):
        alb_lst = []
        for key in retrieval_init['trace gases'].keys():
            if (key == 'CO2'):
                atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']
            if (key == 'CH4'):
                atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']
            if (key == 'H2O'):
                atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']

        m = 0
        while "alb%d" % (m) in retrieval_init['surface'].keys():
            alb_lst.append(x_dic["alb%d" % (m)])
            m = m+1

        x0_lst = np.array(list(x_dic.values()))  # dictionary values -> list -> numpy

        # Calculate surface data

        surface.get_albedo_poly(alb_lst)

        # Calculate nonscattered forward model
        fwd = libRT.nonscat_fwd_model(isrf, retrieval_init['solar irradiance'],
                                      atm,  optics, surface, measurement['mu0'],
                                      measurement['muv'], dev)

        ytilde = measurement['ymeas'] - fwd['rad']  # Difference between forward model and measured spectrum

        # Calculate convoluted Jacobians
        Kmat = np.zeros([len(fwd['rad']), len(x_dic)])

        l = 0
        for value in k_dic.values():
            Kmat[:, l] = fwd[value]
            l += 1

        # Calculated least square solution
        Syinv = np.eye(fwd['rad'].size)*1./np.diag(measurement['Smeas'])
#        Syinv = np.linalg.inv(measurement['Smeas'])           # inverse of covariance matrix of the measurement
        # covariance matrix of the estimated least square solution
        Sx = np.linalg.inv(np.dot(Kmat.T, np.dot(Syinv, Kmat)))
        Gain = np.dot(Sx, np.dot(Kmat.T, Syinv))               # gain matrix
        xstat = np.dot(Gain, ytilde)                           # least square solution


        # print(x0_lst[0:1])
        # print('==========================================')
        # print(iteration, xstat)
        
        # fig = plt.figure(figsize=(10, 8), dpi=100)
        # plt.plot(measurement['ymeas'], color='blue', label='l1b')
        # plt.plot(fwd['rad'], color='green', label='fwd')
        # plt.xlabel('spec index')
        # plt.ylabel('radiance')
        # plt.legend()        
        # sys.exit()

        x_lst_precision = []
        for m in range(0, len(xstat)):
            x_lst_precision.append(np.sqrt(Sx[m, m]))

        k = 0
        x_lst = []
        for key in x_dic.keys():
            # remember that we estimated x-x0, calculate x
            x_dic[key] = x0_lst[k]+xstat[k]
            x_lst.append(x0_lst[k]+xstat[k])
            g_dic[key] = Gain[k, :]
            k += 1

        Chi_sqrt.append(ytilde.T@Syinv @ ytilde/(len(fwd["rad"])-len(x0_lst)))

        # Define convergence criteria or if iterations are too large
        if iteration > 2:
            if (np.abs(Chi_sqrt[iteration]-Chi_sqrt[iteration-1]) < retrieval_init['chi2 limit']):
                convergence = True
                break
        iteration = iteration+1

    output = {}
    output['chi2'] = Chi_sqrt[iteration-1]
    output['convergence'] = convergence
    output['number_iter'] = iteration
    # define output product, first update all parameter
    for key in retrieval_init['trace gases'].keys():
        if (key == 'CO2'):
            atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']
        if (key == 'CH4'):
            atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']
        if (key == 'H2O'):
            atm.__getattribute__(key)[:] = x_dic[key]*retrieval_init['trace gases'][key]['ref_profile']

    # calculate column mixing ratio, precision and albedo values
    for l, key in enumerate(retrieval_init['trace gases'].keys()):
        scaling = retrieval_init['trace gases'][key]['scaling']
        output['X'+key] = sum(atm.__getattribute__(key))/sum(atm.air)
        output['X'+key+' precision'] = output['X'+key] * x_lst_precision[l]           # ppm
        ref_mixing_ratio = np.sum(retrieval_init['trace gases'][key]['ref_profile'])/sum(atm.air)/scaling
        output['gain_X'+key] = g_dic[key] * ref_mixing_ratio

    m = 0
    while 'alb%d' % (m) in retrieval_init['surface'].keys():
        output['alb%d' % (m)] = x_dic["alb%d" % (m)]
        m = m+1

    # Calculate column averaging kernels

    nlay = atm.zlay.size
    output['XCO2 col avg kernel'] = np.zeros(nlay)
    output['XCH4 col avg kernel'] = np.zeros(nlay)
    output['XH2O col avg kernel'] = np.zeros(nlay)

    # we have to use k_dic because of the labeling of the forward model derivatives.
    # Something to change in a next version. I get made about the differnet labeling

    for ispec, value in enumerate(k_dic.values()):
        if(value == 'molec_07'):  # CO2
            col = np.sum(retrieval_init['trace gases']['CO2']['ref_profile'])
            for klay in range(nlay):
                delta_prof = retrieval_init['trace gases']['CO2']['ref_profile'][klay]/col
                output['XCO2 col avg kernel'][klay] = np.sum(Gain[ispec, :]*fwd['layer_'+value][:, klay])*delta_prof
        if(value == 'molec_32'):  # CH4
            col = np.sum(retrieval_init['trace gases']['CH4']['ref_profile'])
            for klay in range(nlay):
                delta_prof = retrieval_init['trace gases']['CH4']['ref_profile'][klay]/col
                output['XCH4 col avg kernel'][klay] = np.sum(Gain[ispec, :]*fwd['layer_'+value][:, klay])*delta_prof
        if(value == 'molec_01'):  # H2O
            col = np.sum(retrieval_init['trace gases']['H2O']['ref_profile'])
            for klay in range(nlay):
                delta_prof = retrieval_init['trace gases']['H2O']['ref_profile'][klay]/col
                output['XH2O col avg kernel'][klay] = np.sum(Gain[ispec, :]*fwd['layer_'+value][:, klay])*delta_prof

    return output
