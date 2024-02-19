import numpy as np
import os
import sys
import inspect
import logging
import yaml

import netCDF4 as nc

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from lib import libDOAS


def ifdoe_run(cfg):
    '''
     Routine which contains the IFDOE processing chain, i.e
     the main sections.
    '''

    verboseLevel = 4

    # only one set of ref spec for one pixel
    groundPixel = 0


    logger = logging.getLogger()
    
    logger.info('Start IFDOE processing')
    
    cfg = libDOAS.AddConfig(cfg)

    # A.2  Assemble list of parameters
    
    parameterNames = libDOAS.IFDOEparameters(cfg)


    # A.3  Read ref.spec. files & rescale wavelengths to [-1:+1] 

    logging.info('Reading ref.spec.')

    # read in ascii ref.spec. for given ground pixel

    IFDOERefSpec = libDOAS.ReadRefSpecAsc(cfg,groundPixel)

    # rescale wvl to [-1:+1] for fit window
    IFDOERefSpec['solarWvlScaled'] = libDOAS.ScaleWvl(IFDOERefSpec['wvl'],cfg['fit_window'])
    IFDOERefSpec['xsWvlScaled'] = IFDOERefSpec['solarWvlScaled'] 

    
    # A.4  Setup loop bounds

    scanN,pxlN,spectralN = libDOAS.getDimensionsRad(cfg['input']['rad'])

    if 'alt' in cfg:
        scanBeg = cfg['alt']['start']
        scanEnd = cfg['alt']['stop']
    else:
        scanBeg = 0
        scanEnd = scanN - 1

    if 'act' in cfg:
        pxlBeg = cfg['act']['start']
        pxlEnd = pxlRange['stop']
    else:
        pxlBeg = 0
        pxlEnd = pxlN - 1


    # A.5  Initialise output

    results = {}
    results['errorFlag'] = np.full((scanN,pxlN),np.nan)

    # B)  Solar spectrum
    # ------------------


    # B.1  Read spectrum & error

    irrWvl, irr, irrError, irrFlag = libDOAS.ReadIrrSGM(cfg['input']['irr'], pxlN)

    # B.2  Filter for irradiance flags
    if (irrFlag!=0).any():
        logging.info('Irr spectra: groundpixel(s) {} flagged'.format(np.where((irrFlag!=0).any(axis=1))[0][:]))
        results['errorFlag'][:,(irrFlag!=0).any(axis=1)] = 1

    # B.3  Calibrate spectrum

    if cfg['wlcal_solar']:
        logging.info('=== Calibrating solar spectrum')

        # calibrate spectrum for each ground pixel
        irrWvlCal       = np.full((pxlN,irrWvl.shape[-1]),np.nan)
        irrCalWs        = np.full((pxlN),np.nan)
        irrCalWsSigma   = np.full((pxlN),np.nan)
        irrCalWq        = np.full((pxlN),np.nan)
        irrCalWqSigma   = np.full((pxlN),np.nan)
        irrCalChiSquare = np.full((pxlN),np.nan)
        irrCalConverged = np.full((pxlN),np.nan)

        for ipxl in range(pxlBeg,pxlEnd+1,1):
            irrWvlCal[ipxl,:], irrCalWs[ipxl], irrCalWsSigma[ipxl], irrCalWq[ipxl], irrCalWqSigma[ipxl], irrCalChiSquare[ipxl], irrCalConverged[ipxl]= libDOAS.irrCal(irr[ipxl,:],irrError[ipxl,:],irrWvl[ipxl,:],cfg,IFDOERefSpec,groundPixel,verboseLevel)
            breakpoint()
            try:
                irrWvlCal[ipxl,:], irrCalWs[ipxl], irrCalWsSigma[ipxl], irrCalWq[ipxl], irrCalWqSigma[ipxl], irrCalChiSquare[ipxl], irrCalConverged[ipxl]= libDOAS.irrCal(irr[ipxl,:],irrError[ipxl,:],irrWvl[ipxl,:],cfg,IFDOERefSpec,groundPixel,verboseLevel)
            except:
                logging.error(f'Pixel {ipxl} not converged')
                irrCalConverged[ipxl] = 0

        # flag non-converged spectra
        if (irrCalConverged==0).any():
            logging.info('Irr wvl cal: groundpixel(s) {} not converged'.format(np.where(irrCalConverged==0)))
            results['errorFlag'][:,irrCalConverged==0] = 1

        else:
            irrWvlCal       = np.full((pxlN,irrWvl.shape[-1]),np.nan)
            irrCalWs = irrCalWsSigma = irrCalWq = irrCalWqSigma = irrCalChiSquare = irrCalConverged = np.full((pxlN),np.nan)

        # B.4  Write calibration results to output_file

        irrCalResults = ['irrWvlCal','irrCalWs','irrCalWsSigma','irrCalWq','irrCalWqSigma','irrCalChiSquare']
        for parameter in irrCalResults:
            results[parameter] = eval(parameter)

    # B.5  Rescale wavelengths to [-1:+1]

    if cfg['wlcal_solar']:
        irrWvlScaled = libDOAS.ScaleWvl(irrWvlCal,cfg['fit_window'])
    else:
        irrWvlScaled = libDOAS.ScaleWvl(irrWvl,cfg['fit_window'])

    breakpoint()
    
  # C)  Loop over scanlines
  # -----------------------

  # C.1  Set up arrays for the scanline results

    resultScalarList = ['RMS','DOF','chiSquare','nIter','nOutliers','converged','nValid',
                        'radCalWs','radCalWsSigma','radCalWq','radCalWqSigma', 'radCalChiSquare']


    for name in resultScalarList:
        results[name]           = np.full((scanN,pxlN),np.nan)

    # resultVectorList = ['RModel','RMeas','RRes','RTOA','wvl']

    # for name in resultVectorList:
    #     results[name]           = np.zeros((scanN,pxlN,wvlN))

    for name in parameterNames:
        results[name]           = np.full((scanN,pxlN),np.nan)
        results[name+'Sigma']   = np.full((scanN,pxlN),np.nan)

    if (ProgKind == 'trop') or (ProgKind == 'co2m'):
        results['amfGeo'] = np.full((scanN,pxlN),np.nan)
        results['NO2Geo'] = np.full((scanN,pxlN),np.nan)
        results['NO2GeoSigma'] = np.full((scanN,pxlN),np.nan)
        results['R_NO2'] = np.full((scanN,pxlN),np.nan)
        results['R_NO2_precision'] = np.full((scanN,pxlN),np.nan)
        results['R_O2O2'] = np.full((scanN,pxlN),np.nan)
        results['R_O2O2_precision'] = np.full((scanN,pxlN),np.nan)

    if exportSpectra:  
        results['wvl'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_meas'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_model'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_res'] = np.full((scanN,pxlN,spectralN),np.nan)


    for iscan in range(scanBeg,scanEnd+1,1):
        startTime = datetime.now()

        # C.2  Read earth spectra & error

        if ProgKind == 'trop':
            radWvlScan, radScan, radErrorScan, radFlagScan, szaScan, vzaScan = IFDOE_func.ReadRadTrop(radFile,iscan)

        elif ProgKind == 'co2m':
            radWvlScan, radScan, radErrorScan, radFlagScan, szaScan, vzaScan = IFDOE_func.ReadRadCO2M(radFile,iscan,rng)

        elif ProgKind == 'asc' and irrFile != None:
            radWvl, rad, radError = IFDOE_func.ReadSpecAscii(radFile)


        if (ProgKind == 'trop') or (ProgKind == 'co2m'):

            # Filter for radiance flags
            if (radFlagScan>2).any():

                # allow eclipse and sun glint flag
                IFDOE_func.printMessage(verboseLevel,'Rad spectrum: scanline {}, groundpixel(s) {} flagged'.format(iscan,np.where(radFlagScan>2)[0][:]))
                results['errorFlag'][iscan,radFlagScan>2] = 1


        
        # C.3  Loop over ground pixels along a scanline
        
        for ipxl in range(pxlBeg,pxlEnd+1,1):

            if irrFile != None:

                # start = time.time()

              # D)  Earth spectra
              # -----------------

                # D.1  Get spectrum & error

                if (ProgKind == 'trop') or (ProgKind == 'co2m'):

                    rad = radScan[ipxl,:]
                    radError = radErrorScan[ipxl,:]
                    radWvl = radWvlScan[ipxl]
                    sza = szaScan[ipxl]
                    vza = vzaScan[ipxl]


                    # skip pixel if sza > sza_limit
                    if sza > cfg['sza_limit']:
                        results['errorFlag'][iscan,ipxl] = 1
                        IFDOE_func.printMessage(verboseLevel,'Scanline {}, ground pixel {}: SZA:{} > SZA_limit {}'.format(iscan,ipxl,sza,cfg['sza_limit']))
                        continue

                    #  skip if error flag is already raised
                    if results['errorFlag'][iscan,ipxl] == 1 :
                        IFDOE_func.printMessage(verboseLevel,'Skipping scanline {}, pixel {}: error flag'.format(iscan,ipxl))
                        continue


                # D.2  Calibrate spectrum
                if cfg['wlcal_earth']:
                    IFDOE_func.printMessage(verboseLevel,'=== Calibrating Earth spectrum')

                    try:
                        if cfg['ref_spectra_source'] == 'gaussian':
                            radWvlCal, radCalWs, radCalWsSigma, radCalWq, radCalWqSigma, radCalChiSquare, radCalConverged = IFDOE_cal.radCal(rad,radError,radWvl,cfg,IFDOERefSpec,0,verboseLevel)

                        else:
                            radWvlCal, radCalWs, radCalWsSigma, radCalWq, radCalWqSigma, radCalChiSquare, radCalConverged = IFDOE_cal.radCal(rad,radError,radWvl,cfg,IFDOERefSpec,ipxl,verboseLevel)
                    except:
                        radCalConverged = 0

                    # flag non-converged spectra
                    if radCalConverged!=1:
                        IFDOE_func.printMessage(verboseLevel,'Rad wvl cal: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                        results['errorFlag'][iscan,ipxl] = 1
                        continue

                else:
                    radWvlCal       = np.full((irrWvl.shape[-1]),np.nan)
                    radCalWs = radCalWsSigma = radCalWq = radCalWqSigma = radCalChiSquare = radCalConverged = np.nan


                # D.3  Rescale wavelengths to [-1:+1]

                if cfg['wlcal_earth']:
                    radWvlScaled = IFDOE_func.ScaleWvl(radWvlCal,cfg['fit_window'])
                else:
                    radWvlScaled = IFDOE_func.ScaleWvl(radWvl,cfg['fit_window'])



              # E)  Reflectance
              # ---------------
            
                # # E.1  Interpolate earth grid to solar grid

                # # Interpolate Earth to solar using spline
                # radSolarGrid = SplineInterp1D(radWvlScaled,rad,commonWvl)
                # radErrorSolarGrid = SplineInterp1D(radWvlScaled,radError,commonWvl)

                # splineplotx = np.linspace(-1.2,2.1,10000)
                # splineploty = SplineInterp1D(radWvlScaled,rad,splineplotx)


                # import matplotlib.pyplot as plt
                # plt.figure()
                # plt.plot(radWvlScaled, rad, 'x', label='original')
                # plt.plot(irrWvlScaled, radSolarGrid, 'x', label='regridded')
                # plt.plot(splineplotx, splineploty, '-', label='spline')
                # plt.legend()
                # plt.show()


                # # E.1  Interpolate solar grid to earth grid using ref_solar

                irrRegridded, irrErrorRegridded = IFDOE_func.SolarToEarthWvl(IFDOERefSpec,radWvlScaled,irrWvlScaled,irr,irrError,ipxl,ProgKind)

                commonWvl = radWvlScaled


                # import matplotlib.pyplot as plt
                # plt.figure()
                # plt.plot(irrWvlScaled[ipxl,:], irr[ipxl,:], 'x', label='original')
                # plt.plot(radWvlScaled, irrRegridded, 'x', label='regridded')
                # plt.plot(IFDOERefSpec['solarWvlScaled'],IFDOERefSpec['solar'][0,:], label='ref_solar')
                # plt.legend()
                # plt.title('irradiance')
                # plt.show()


                # E.2  Compute reflectance & error (do NOT add pi/mu0 term)

                refl = rad/irrRegridded
                reflError = 1/irrRegridded * np.sqrt( radError**2 + (irrErrorRegridded**2 * refl**2) )


                # test with pi/mu0
                # refl = radSolarGrid/irr * np.pi/np.cos(25.17924*np.pi/180.0)


                # OR: Read reflectance & error (may have pi/mu0 term)
            # When no irriadiance file is present, assume rad file is refl file
            if ProgKind == 'asc'  and irrFile == None:

                reflWvl, refl, reflError = IFDOE_func.ReadSpecAscii(radFile)

                # rescale wavelengths to [-1:+1]
                # set refl grid as common grid
                commonWvl = IFDOE_func.ScaleWvl(reflWvl,cfg['fit_window'])


            # E.3  Put ref.spec. at common grid using spline; cf. B.6

            if cfg['ref_spectra_source'] == 'gaussian':
                IFDOERefSpecRegrid = IFDOE_func.InterpRefSpecGrid(IFDOERefSpec,commonWvl,ProgKind,cfg,0)
            else:
                IFDOERefSpecRegrid = IFDOE_func.InterpRefSpecGrid(IFDOERefSpec,commonWvl,ProgKind,cfg,ipxl)

    
          # F)  DOAS fit
          # ------------

            IFDOE_func.printMessage(verboseLevel,'=== DOAS fit')

            # F.1  Setup the model details

            # clip spectra to fit window 

            commonWvlFit = commonWvl[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]
            reflFit = refl[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]
            reflErrorFit = reflError[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]


            # plt.figure()
            # plt.plot(commonWvlFit, reflFit)
            # plt.title('reflectance')
            # plt.show()



            IFDOERefSpecFit = {}
            for RefSpec in IFDOERefSpecRegrid:
                IFDOERefSpecFit[RefSpec] = IFDOERefSpecRegrid[RefSpec][(commonWvl>=-1.0 ) & (commonWvl<=1.0)]

            # setup model

            prior = np.asarray([cfg['prior'][key] for key in parameterNames])
            priorCov = np.diag(np.asarray([cfg['prior_error'][key] for key in parameterNames])**2)

            try:
                model = IFDOE_model.IFDOEmodel(prior=prior, priorCov=priorCov, 
                  otherModelParam=None, 
                  parameterNames=parameterNames,
                  verbose=verboseLevel-1,
                  observation=reflFit, 
                  observationError=reflErrorFit,
                  independentVariable=commonWvlFit,
                  cfg=cfg,
                  IFDOErefspec=IFDOERefSpecFit)

                
                # F.2  Perform fit

                oe = IFDOE_OE.OE(model=model, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
                
                # pdb.set_trace()
                oe()

            except:
                results['errorFlag'][iscan,ipxl] = 1
                IFDOE_func.printMessage(verboseLevel,'DOAS : scanline {}, groundpixel {} unexpected error'.format(iscan,ipxl))

            IFDOE_func.printMessage(verboseLevel,( 'converged >',oe.answer['converged'] ))

            # flag non-converged fit
            if oe.answer['converged']==False:
                results['errorFlag'][iscan,ipxl] = 1
                IFDOE_func.printMessage(verboseLevel,'DOAS: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                continue



            if cfg['spike_removal']:
            # F.3  Spike removal -> re-do E.2, but only once

                # detect outliers
                
                wvl = IFDOE_func.ScaleWvl(commonWvlFit,cfg['fit_window'])
                RRes = model.observation - model.modelCalculation
                Q1 = np.quantile(RRes, 0.25)
                Q3 = np.quantile(RRes, 0.75)
                Qf = 3.0
                outerFenceUpper = Q3+ Qf*(Q3-Q1)
                outerFenceLower = Q1- Qf*(Q3-Q1)

                outlierIndex = (RRes>outerFenceUpper) | (RRes<outerFenceLower)
                outlierWvl = wvl[outlierIndex]
                outlierRRes = RRes[outlierIndex]
                nOutliers = outlierIndex.sum()

                outliers = False
                if nOutliers > 0 :
                    outliers = True
                if outliers:
                    IFDOE_func.printMessage(verboseLevel,('{} outliers detected at wvl: '.format(outlierIndex.sum()),outlierWvl))
                else:
                    IFDOE_func.printMessage(verboseLevel,('No outliers detected.'))

                if nOutliers > cfg['max_outliers']:
                    IFDOE_func.printMessage(verboseLevel,('Scanline {}, ground pixel {} : too many outliers ({}), raising error flag'.format(iscan,ipxl,nOutliers)))
                    results['errorFlag'][iscan,ipxl] = 1
                    continue


                if outliers:
                    # remove outliers
                    reflFit = reflFit[~outlierIndex]
                    reflErrorFit = reflErrorFit[~outlierIndex]
                    commonWvlFit = commonWvlFit[~outlierIndex]

                    for RefSpec in IFDOERefSpecFit:
                        IFDOERefSpecFit[RefSpec] = IFDOERefSpecFit[RefSpec][~outlierIndex]

                    # Perform DOAS fit once again
                    try:

                        model = IFDOE_model.IFDOEmodel(prior=prior, priorCov=priorCov, 
                          otherModelParam=None, 
                          parameterNames=parameterNames,
                          verbose=verboseLevel-1,
                          observation=reflFit, 
                          observationError=reflErrorFit,
                          independentVariable=commonWvlFit,
                          cfg=cfg,
                          IFDOErefspec=IFDOERefSpecFit)

                        oe = IFDOE_OE.OE(model=model, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
                        oe()

                    except:
                        results['errorFlag'][iscan,ipxl] = 1
                        IFDOE_func.printMessage(verboseLevel,'DOAS after outlier removal: scanline {}, groundpixel {} unexpected error'.format(iscan,ipxl))

                    # flag non-converged fit
                    if oe.answer['converged']==False:
                        results['errorFlag'][iscan,ipxl] = 1
                        IFDOE_func.printMessage(verboseLevel,'DOAS after outlier removal: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                        continue
            else:
                nOutliers = np.nan


        # F.4  Determine diagnostics: RMS error, degrees of freedom, etc.

            RMS = model.rmse
            DOF = oe.answer['degrees of freedom']
            chiSquare = oe.costFunction
            nIter = oe.answer['number of iterations']
            converged = oe.answer['converged']
            nValid = commonWvlFit.shape[0]

            wvl = IFDOE_func.DeScaleWvl(commonWvlFit,cfg['fit_window'])
            RModel = model.modelCalculation
            RMeas = model.observation
            RRes = RMeas - RModel


            # scale precision with chisquare reduced
            chiSquare_reduced = chiSquare/(nValid - DOF)
            precision = np.diag(np.sqrt(model.covariance * chiSquare_reduced))

            if ProgKind=='asc':

                polyCoeff = np.zeros((cfg['polynomial_coefs']))
                polyCoeffSigma = np.zeros((cfg['polynomial_coefs']))

                for i in range(len(polyCoeff)):
                    polyCoeff[i] = model.stateVector[i]
                    polyCoeffSigma[i] = precision[i]

                intensityCoeff = np.zeros((cfg['intensity_coefs']))
                intensityCoeffSigma = np.zeros((cfg['intensity_coefs']))

                for i in range(len(intensityCoeff)):
                    j= i+cfg['polynomial_coefs']
                    intensityCoeff[i] = model.stateVector[j]
                    intensityCoeffSigma[i] = precision[j]

                gasColumns = {}
                gasColumnsSigma = {}

                for i,gas in enumerate(cfg['trace_gas_list']):
                    i+= cfg['polynomial_coefs']+cfg['intensity_coefs']
                    gasColumns[gas] = model.stateVector[i] 
                    gasColumnsSigma[gas] = precision[i]

                if cfg['ring_fit_term'] == 'Iring':
                    i = cfg['polynomial_coefs']+cfg['intensity_coefs']+cfg['nr_trace_gases']
                    Cring = model.stateVector[i] 
                    CringSigma = precision[i]

                elif cfg['ring_fit_term'] == 'Dring':
                    i = cfg['polynomial_coefs']+cfg['intensity_coefs']+cfg['nr_trace_gases']
                    Dring = model.stateVector[i] 
                    DringSigma = precision[i]

            # F.5  Determine as for TROPOMI: R_toa = P * (1+C_ring)
                # Without pi/mu0

                RTOA = np.zeros((len(wvl)))
                for i in range(len(polyCoeff)):
                    RTOA += polyCoeff[i] * commonWvlFit**i

                if cfg['ring_fit_term'] == 'Iring':
                    RTOA *= (1 + Cring)

                # print(RTOA[np.abs(wvl-440).argmin()]*(np.pi/np.cos(np.pi*25.179/180)))

        # F.6  Show/plot results

            # Print the results
            # -----------------

            if detailLevel == 2:

                # The "standard" overview or results (and model)   
                print(oe)
                print(model)

                # # All input/output connected to 'oe'

                print( '' )
                print( '*************' )
                print( '' )

                al = oe.answer.keys()
                for a in (al) :
                    print( a,' > ',oe.answer[a] )

            if verboseLevel > 1 :
                print('\n========================================================')
                print('Summary fitted parameters')
                print('========================================================')
                print('{:<10}{:>15}{:>15}'.format('(parameter)','(value)','(std)'))
                for i in range(len(model.stateVector)):
                    print('{:<10s}{:>15.6E}{:>15.6E}'.format(model.parameterNames[i], model.stateVector[i], precision[i]))

            if verboseLevel == 3 :

                # Plot Rmod, Rmeas and residual
                # -------------------

                plt.figure()
                plt.plot(wvl, RModel, '-', label='Model')
                plt.plot(wvl, RMeas, '--', label='Measurement')
                plt.xlabel('wavelength [nm]')
                plt.ylabel('reflectance [-]')
                plt.legend()
                plt.draw()

                plt.figure()
                plt.plot(wvl, RRes, label='Residual')
                plt.xlabel('wavelength [nm]')
                plt.ylabel('residual [-]')
                plt.legend()
                plt.draw()

                plt.show(block=False)

                s = input('\n>>> Hit enter to stop : ')
                plt.close("all")

        # C.4  Write scanline results & debug data to output

            # scalar
            for parameter in resultScalarList:
                results[parameter][iscan,ipxl] = eval(parameter)

            # # vector
            # for parameter in resultVectorList:
            #     results[parameter][iscan,ipxl,:] = eval(parameter)

            # if geometry info available, scale poly coeff
            if (ProgKind == 'trop') or (ProgKind == 'co2m'):
                szaFactor = np.pi / np.cos(sza*np.pi/180.0)
            else:
                szaFactor = 1.0

            # fit parameters
            for i,parameter in enumerate(parameterNames):

                # scale poly coeff with pi/mu0
                if 'P' in parameter and (ProgKind=='trop' or ProgKind=='co2m'):
                    results[parameter][iscan,ipxl] = model.stateVector[i] * szaFactor
                    results[parameter+'Sigma'][iscan,ipxl] = precision[i] * szaFactor
                else:
                    results[parameter][iscan,ipxl] = model.stateVector[i] 
                    results[parameter+'Sigma'][iscan,ipxl] = precision[i]

            # if geometry info available calculate amf geo and geo column
            if (ProgKind == 'trop') or (ProgKind == 'co2m'):
                amfGeo = 1/np.cos(sza*np.pi/180.0) + 1/np.cos(vza*np.pi/180.0) 
                results['amfGeo'][iscan,ipxl] = amfGeo
                results['NO2Geo'][iscan,ipxl] = results['NO2'][iscan,ipxl] / amfGeo
                results['NO2GeoSigma'][iscan,ipxl] = results['NO2Sigma'][iscan,ipxl] / amfGeo

            

            if (ProgKind == 'trop') or (ProgKind == 'co2m'):

                # continuum reflectance at 440 nm for NO2:

                if cfg['fit_window'][1] < 477.0:

                    iwvl = np.abs(wvl-440.0).argmin()
                    results['R_NO2'][iscan,ipxl] = 0.0
                    results['R_NO2_precision'][iscan,ipxl] = 0.0
                    
                    for i in range(cfg['polynomial_coefs']):

                        results['R_NO2'][iscan,ipxl] += model.stateVector[i] *commonWvlFit[iwvl] **i
                        results['R_NO2_precision'][iscan,ipxl] += precision[i] *commonWvlFit[iwvl] **i

                    if cfg['ring_fit_term'] == 'Iring':
                        i = cfg['polynomial_coefs']+cfg['intensity_coefs']+cfg['nr_trace_gases']
                        Cring = model.stateVector[i] 
                        results['R_NO2'][iscan,ipxl] *= (1 + Cring)

                    results['R_NO2'][iscan,ipxl] *= szaFactor
                    results['R_NO2_precision'][iscan,ipxl] *= szaFactor

                # continuum reflectance at 477 nm of O2O2:

                elif cfg['fit_window'][1] > 477.0:

                    iwvl = np.abs(wvl-477.0).argmin()
                    results['R_O2O2'][iscan,ipxl] = 0.0
                    results['R_O2O2_precision'][iscan,ipxl] = 0.0
                    
                    for i in range(cfg['polynomial_coefs']):

                        results['R_O2O2'][iscan,ipxl] += model.stateVector[i] *commonWvlFit[iwvl] **i
                        results['R_O2O2_precision'][iscan,ipxl] += precision[i] *commonWvlFit[iwvl] **i


                    # do not add Cring for O2O2

                    # if cfg['ring_fit_term'] == 'Iring':
                    #     i = cfg['polynomial_coefs']+cfg['intensity_coefs']+cfg['nr_trace_gases']
                    #     Cring = model.stateVector[i] 
                    #     results['R_O2O2'][iscan,ipxl] *= (1 + Cring)

                    results['R_O2O2'][iscan,ipxl] *= szaFactor
                    results['R_O2O2_precision'][iscan,ipxl] *= szaFactor


            if exportSpectra:

                results['wvl'][iscan,ipxl,:len(wvl)] = wvl
                results['R_meas'][iscan,ipxl,:len(wvl)] = RMeas*szaFactor
                results['R_model'][iscan,ipxl,:len(wvl)] = RModel*szaFactor
                results['R_res'][iscan,ipxl,:len(wvl)] = RRes*szaFactor
                    

            # print(time.time()-start)
            results['errorFlag'][iscan,ipxl] = 0


        if ((ProgKind == 'trop') or (ProgKind == 'co2m')) and verboseLevel > 0:
            print('Processed scanline {}/{} in {}s'.format(iscan-scanBeg+1,scanEnd-scanBeg+1,np.round((datetime.now() - startTime).total_seconds(),2) ))


    # H)  Finishing up

    if ProgKind == 'trop':

        # check if output file exists:
        if os.path.exists(outFile) == False:
            # write results to netcdf file (copy from Tropomi NO2 L2 format)
            IFDOE_func.initOutput(outFileFormat, outFile)

            # copy fields from L1B file to output file
            IFDOE_func.copyOutputTrop(radFile,outFile)    

        # write results to output file
        IFDOE_func.writeOutputTrop(cfg,parameterNames,results,outFile,exportSpectra)

    elif ProgKind == 'co2m':

        # check if output file exists:
        if os.path.exists(outFile) == False:

            # write results to netcdf file (copy from CO2M NO2 L2 format)
            IFDOE_func.initOutputXml(outFileFormat, outFile)

            # copy fields from L1B file to output file
            IFDOE_func.copyOutputCO2M(radFile,outFile)    

        # write results to output file
        IFDOE_func.writeOutputCO2M(cfg,parameterNames,results,outFile,exportSpectra)

    if ProgKind != 'asc' and verboseLevel > 0:
            print(f'Output witten to {outFile}')
    

    return


if __name__ == '__main__':
    

    # call with:
    # python l1bl2_no2.py l1bl2_no2.yaml

    # or with logging to file:
    # python l1bl2_no2.py l1bl2_no2.yaml l1bl2_no2.log

    # setup the logging to screen and to file
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

    if len(sys.argv) > 2:
        fh = logging.FileHandler(sys.argv[2], mode='w')
        fh.setLevel(logging.ERROR)
        fh.setFormatter(formatter)

        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter) 

        logging.basicConfig(level=logging.INFO, handlers = [ch,fh])

        logging.info(f'Logging to file: {sys.argv[2]}')

    else:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter) 
        logging.basicConfig(level=logging.INFO,handlers = [ch])

    logging.info(f'Reading config file: {sys.argv[1]}')



    cfg = yaml.safe_load(open(sys.argv[1]))


    ifdoe_run(cfg['doas'])
