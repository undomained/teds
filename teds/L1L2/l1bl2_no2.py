import numpy as np
import os
import sys
import inspect
import logging
import yaml
import netCDF4 as nc
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
import multiprocessing.pool
import time
from scipy.ndimage import gaussian_filter1d
import datetime

from teds.lib import libDOAS, libAMF
from teds.lib.libWrite import writevariablefromname

def conv_irr(sgm_rad_file, fwhm):
    # convolve irradiance with Gaussian ISRF
    with nc.Dataset(sgm_rad_file) as f:

        irr = f['solar_irradiance'][:] # [spectral_bins] - "photons / (nm m2 s)"
        wvl = f['wavelength'][:] # [spectral_bins] - nm

    stepsize = wvl[1]-wvl[0]
    fwhm_step = fwhm/stepsize
    sigma = fwhm_step /np.sqrt(8*np.log(2))
    convolved_irr = gaussian_filter1d(irr, sigma)

    file_out = sgm_rad_file.replace('.nc','_conv_irr.nc')
    # open file
    with nc.Dataset(file_out, mode='w') as output_conv_irr:
        output_conv_irr.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
        output_conv_irr.comment = f'Convolved irradiance with Gaussian FWHM {fwhm} nm'
        output_conv_irr.createDimension('bins_spectral', len(irr))     # spectral axis
        # wavelength
        _ = writevariablefromname(output_conv_irr, 'wavelength', ('bins_spectral',), wvl)
        # solar irradiance
        _ = writevariablefromname(output_conv_irr, 'solarirradiance', ('bins_spectral',), convolved_irr)

    return file_out


def ifdoe_run(logger, cfg):
    '''
     Routine which contains the IFDOE processing chain, i.e
     the main sections.
    '''

    if cfg['debug']['plot']:
        verboseLevel = 3
    elif cfg['debug']['log']:
        verboseLevel = 2
    else:
        verboseLevel = 0

    # only one set of ref spec for one pixel TODO: change when ISRF is known
    groundPixel = 0

    startTime = time.time()
    
    logger.info('=== Start IFDOE processing')
    
    cfg = libDOAS.AddConfig(cfg)

    # A.2  Assemble list of parameters
    
    parameterNames = libDOAS.IFDOEparameters(cfg)


    # A.3  Read ref.spec. files & rescale wavelengths to [-1:+1] 

    logger.info('=== Reading ref.spec.')

    # read in ascii ref.spec. for given ground pixel

    IFDOERefSpec = libDOAS.ReadRefSpecAsc(cfg,groundPixel)

    # rescale wvl to [-1:+1] for fit window
    IFDOERefSpec['solarWvlScaled'] = libDOAS.ScaleWvl(IFDOERefSpec['wvl'],cfg['fit_window'])
    IFDOERefSpec['xsWvlScaled'] = IFDOERefSpec['solarWvlScaled'] 

    
    # A.4  Setup loop bounds

    # TODO rad_file did not exist in full_config.yaml. SHOULD it be cfg['io']['sgm_rad'] ?????
    scanN,pxlN,spectralN = libDOAS.getDimensionsRad(cfg['rad_file'])

    if 'alt' in cfg:
        scanBeg = cfg['alt']['start']
        scanEnd = cfg['alt']['stop']
    else:
        scanBeg = 0
        scanEnd = scanN - 1

    if 'act' in cfg:
        pxlBeg = cfg['act']['start']
        pxlEnd = cfg['act']['stop']
    else:
        pxlBeg = 0
        pxlEnd = pxlN - 1


    # A.5  Initialise output

    results = {}
    results['errorFlag'] = np.full((scanN,pxlN),np.nan)


    # read geometry

    geo = libDOAS.readGeometry(cfg['io']['gm'])

    # B)  Solar spectrum
    # ------------------


    # B.1  Read spectrum & error
    if cfg['irr_from_sgm']:
        # TODO: irr_file does not exist in full_config.yaml. what should this be?
        irrWvl, irr, irrError, irrFlag = libDOAS.ReadIrrSGM(cfg['irr_file'], pxlN)
    else:
        # TODO: not implemented yet
        sys.exit()
        # irrWvl, irr, irrError, irrFlag = libDOAS.ReadIrr(cfg['irr_file'], pxlN)

    # B.2  Filter for irradiance flags
    if (irrFlag!=0).any():
        logger.info('Irr spectra: groundpixel(s) {} flagged'.format(np.where((irrFlag!=0).any(axis=1))[0][:]))
        results['errorFlag'][:,(irrFlag!=0).any(axis=1)] = 1

    # B.3  Calibrate spectrum

    if cfg['wlcal_solar']:
        logger.info('=== Calibrating solar spectrum')

        # calibrate spectrum for each ground pixel
        irrWvlCal       = np.full((pxlN,irrWvl.shape[-1]),np.nan)
        irrCalWs        = np.full((pxlN),np.nan)
        irrCalWsSigma   = np.full((pxlN),np.nan)
        irrCalWq        = np.full((pxlN),np.nan)
        irrCalWqSigma   = np.full((pxlN),np.nan)
        irrCalChiSquare = np.full((pxlN),np.nan)
        irrCalConverged = np.full((pxlN),np.nan)

        for ipxl in range(pxlBeg,pxlEnd+1,1):
            try:
                irrWvlCal[ipxl,:], irrCalWs[ipxl], irrCalWsSigma[ipxl], irrCalWq[ipxl], irrCalWqSigma[ipxl], irrCalChiSquare[ipxl], irrCalConverged[ipxl] = libDOAS.irrCal(irr[ipxl,:],irrError[ipxl,:],irrWvl[ipxl,:],cfg,IFDOERefSpec,groundPixel,verboseLevel)
            except:
                logger.error(f'Pixel {ipxl} not converged')
                irrCalConverged[ipxl] = 0

        # flag non-converged spectra
        if (irrCalConverged==0).any():
            logger.error('Irr wvl cal: groundpixel(s) {} not converged'.format(np.where(irrCalConverged==0)))
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

    
  # C)  Loop over scanlines
  # -----------------------

  # C.1  Set up arrays for the scanline results

    resultScalarList = ['RMS','DOF','chiSquare','nIter','nOutliers','converged','nValid',
                        'radCalWs','radCalWsSigma','radCalWq','radCalWqSigma', 'radCalChiSquare']

    for name in resultScalarList:
        results[name]           = np.full((scanN,pxlN),np.nan)

    for name in parameterNames:
        results[name]           = np.full((scanN,pxlN),np.nan)
        results[name+'Sigma']   = np.full((scanN,pxlN),np.nan)

    results['amfGeo'] = np.full((scanN,pxlN),np.nan)
    results['NO2Geo'] = np.full((scanN,pxlN),np.nan)
    results['NO2GeoSigma'] = np.full((scanN,pxlN),np.nan)
    results['R_NO2'] = np.full((scanN,pxlN),np.nan)
    results['R_NO2_precision'] = np.full((scanN,pxlN),np.nan)
    results['R_O2O2'] = np.full((scanN,pxlN),np.nan)
    results['R_O2O2_precision'] = np.full((scanN,pxlN),np.nan)

    if cfg['export_spectra']:  
        results['wvl'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_meas'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_model'] = np.full((scanN,pxlN,spectralN),np.nan)
        results['R_res'] = np.full((scanN,pxlN,spectralN),np.nan)

    # open file before threading, to ensure thread safety
    ncRad = nc.Dataset(cfg['rad_file'])

    logger.info('=== Starting radiance calibration and DOAS')

    # scanline loop, inside internal function, so only the scanline index and mutable results dict have to be passed
    def process_scanline(iscan,results):

        startTimeScanline = time.time()

        # C.2  Read earth spectra & error

        if cfg['rad_from_sgm']:
            radWvlScan, radScan, radErrorScan, radFlagScan = libDOAS.ReadRadSGM(ncRad, iscan)
        else:
            radWvlScan, radScan, radErrorScan, radFlagScan = libDOAS.ReadRad(ncRad, iscan)

        # C.3  Loop over ground pixels along a scanline

        for ipxl in range(pxlBeg,pxlEnd+1,1):

            # D)  Earth spectra
            # -----------------
            # D.1  Get spectrum & error

            rad = radScan[ipxl,:]
            radError = radErrorScan[ipxl,:]
            radWvl = radWvlScan[ipxl,:]
            
            #  skip if error flag is already raised
            if results['errorFlag'][iscan,ipxl] == 1 :
                logger.error('Skipping scanline {}, pixel {}: error flag'.format(iscan,ipxl))
                continue


            # D.2  Calibrate spectrum
            if cfg['wlcal_earth']:
                logger.debug('=== Calibrating Earth spectrum')
                
                try:
                    radWvlCal, radCalWs, radCalWsSigma, radCalWq, radCalWqSigma, radCalChiSquare, radCalConverged = libDOAS.radCal(rad,radError,radWvl,cfg,IFDOERefSpec,0,verboseLevel)
                except:
                    radCalConverged = 0

                # flag non-converged spectra
                if radCalConverged!=1:
                    logger.error('Rad wvl cal: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                    results['errorFlag'][iscan,ipxl] = 1
                    continue

            else:
                radWvlCal       = np.full((irr.shape[-1]),np.nan)
                radCalWs = radCalWsSigma = radCalWq = radCalWqSigma = radCalChiSquare = radCalConverged = np.nan


            # D.3  Rescale wavelengths to [-1:+1]

            if cfg['wlcal_earth']:
                radWvlScaled = libDOAS.ScaleWvl(radWvlCal,cfg['fit_window'])
            else:
                radWvlScaled = libDOAS.ScaleWvl(radWvl,cfg['fit_window'])

            # E)  Reflectance
            # ---------------
            # # E.1  Interpolate solar grid to earth grid using ref_solar
            # irrRegridded, irrErrorRegridded = libDOAS.SolarToEarthWvl(IFDOERefSpec,radWvlScaled,irrWvlScaled,irr,irrError,ipxl)

            # # E.1  Interpolate solar grid to earth grid
            irrRegridded, irrErrorRegridded = libDOAS.SolarToEarthWvl_simple(radWvlScaled,irrWvlScaled,irr,irrError,ipxl)

            # plt.figure()
            # plt.plot(radWvlScaled, irrRegridded, label='Irr regrid')
            # plt.plot(irrWvlScaled[ipxl,:], irr[ipxl,:], label='Irr')
            # plt.plot(radWvlScaled, rad/0.022, label='Rad')
            # plt.legend()
            # plt.show()

            commonWvl = radWvlScaled

            # E.2  Compute reflectance & error (do NOT add pi/mu0 term)
            refl = rad/irrRegridded
            reflError = 1/irrRegridded * np.sqrt( radError**2 + (irrErrorRegridded**2 * refl**2) )


            # E.3  Put ref.spec. at common grid using spline; cf. B.6
            IFDOERefSpecRegrid = libDOAS.InterpRefSpecGrid(IFDOERefSpec,commonWvl,cfg,groundPixel)

            # F)  DOAS fit
            # ------------
            logger.debug('=== DOAS fit')

            # F.1  Setup the model details

            # clip spectra to fit window 
            commonWvlFit = commonWvl[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]
            reflFit = refl[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]
            reflErrorFit = reflError[(commonWvl>=-1.0 ) & (commonWvl<=1.0)]

            IFDOERefSpecFit = {}
            for RefSpec in IFDOERefSpecRegrid:
                IFDOERefSpecFit[RefSpec] = IFDOERefSpecRegrid[RefSpec][(commonWvl>=-1.0 ) & (commonWvl<=1.0)]

            # setup model
            prior = np.asarray([cfg['prior']['doas'][key][0] for key in parameterNames])
            priorCov = np.diag(np.asarray([cfg['prior']['doas'][key][1] for key in parameterNames])**2)

            model = libDOAS.IFDOEmodel(prior=prior, priorCov=priorCov, 
                otherModelParam=None, 
                parameterNames=parameterNames,
                verbose=verboseLevel-1,
                observation=reflFit, 
                observationError=reflErrorFit,
                independentVariable=commonWvlFit,
                IFDOEconfig=cfg,
                IFDOErefspec=IFDOERefSpecFit)

            
            # F.2  Perform fit
            try:
                oe = libDOAS.OE(model=model, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
                
                oe()
                logger.debug(f"converged > {oe.answer['converged']}")

                # flag non-converged fit
                if oe.answer['converged']==False:
                    results['errorFlag'][iscan,ipxl] = 1
                    logger.error('DOAS: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                    continue

            except:
                results['errorFlag'][iscan,ipxl] = 1
                logger.error('DOAS : scanline {}, groundpixel {} unexpected error'.format(iscan,ipxl))
                continue

            
            if cfg['spike_removal']:
            # F.3  Spike removal -> re-do E.2, but only once
                # detect outliers   
                wvl = libDOAS.ScaleWvl(commonWvlFit,cfg['fit_window'])
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
                    logger.debug('{} outliers detected at wvl: {}'.format(outlierIndex.sum(),outlierWvl))
                else:
                    logger.debug('No outliers detected.')

                if nOutliers > cfg['max_outliers']:
                    logger.error('Scanline {}, ground pixel {} : too many outliers ({}), raising error flag'.format(iscan,ipxl,nOutliers))
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

                        model = libDOAS.IFDOEmodel(prior=prior, priorCov=priorCov, 
                            otherModelParam=None, 
                            parameterNames=parameterNames,
                            verbose=verboseLevel-1,
                            observation=reflFit, 
                            observationError=reflErrorFit,
                            independentVariable=commonWvlFit,
                            IFDOEconfig=cfg,
                            IFDOErefspec=IFDOERefSpecFit)

                        oe = libDOAS.OE(model=model, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
                        oe()

                        # flag non-converged fit
                        if oe.answer['converged']==False:
                            results['errorFlag'][iscan,ipxl] = 1
                            logger.error('DOAS after outlier removal: scanline {}, groundpixel {} not converged'.format(iscan,ipxl))
                        
                    except:
                        results['errorFlag'][iscan,ipxl] = 1
                        logger.error('DOAS after outlier removal: scanline {}, groundpixel {} unexpected error'.format(iscan,ipxl))
            
            else:
                nOutliers = np.nan

        # F.4  Determine diagnostics: RMS error, degrees of freedom, etc.

            RMS = model.rmse
            DOF = oe.answer['degrees of freedom']
            chiSquare = oe.costFunction
            nIter = oe.answer['number of iterations']
            converged = oe.answer['converged']
            nValid = commonWvlFit.shape[0]

            wvl = libDOAS.DeScaleWvl(commonWvlFit,cfg['fit_window'])
            RModel = model.modelCalculation
            RMeas = model.observation
            RRes = RMeas - RModel

            # scale precision with chisquare reduced
            chiSquare_reduced = chiSquare/(nValid - DOF)
            precision = np.diag(np.sqrt(model.covariance * chiSquare_reduced))

        # F.6  Show/plot results

            # Print the results
            # -----------------

            # The "standard" overview or results (and model)   
            logger.debug(oe)
            logger.debug(model)

            # # All input/output connected to 'oe'

            logger.debug( '' )
            logger.debug( '*************' )
            logger.debug( '' )

            al = oe.answer.keys()
            for a in (al) :
                logger.debug(f'{a} > {oe.answer[a]}')


            logger.debug('========================================================')
            logger.debug('Summary fitted parameters')
            logger.debug('========================================================')
            logger.debug('{:<10}{:>15}{:>15}'.format('(parameter)','(value)','(std)'))
            for i in range(len(model.stateVector)):
                logger.debug('{:<10s}{:>15.6E}{:>15.6E}'.format(model.parameterNames[i], model.stateVector[i], precision[i]))

            if cfg['debug']['plot']:

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

            # scale poly coeff
            szaFactor = np.pi / np.cos(geo['sza'][iscan,ipxl]*np.pi/180.0)

            # fit parameters
            for i,parameter in enumerate(parameterNames):
                if 'P' in parameter:
                    # scale poly coeff with pi/mu0
                    results[parameter][iscan,ipxl] = model.stateVector[i] * szaFactor
                    results[parameter+'Sigma'][iscan,ipxl] = precision[i] * szaFactor
                else:
                    results[parameter][iscan,ipxl] = model.stateVector[i] 
                    results[parameter+'Sigma'][iscan,ipxl] = precision[i]

            # if geometry info available calculate amf geo and geo column
            amfGeo = 1/np.cos(geo['sza'][iscan,ipxl]*np.pi/180.0) + 1/np.cos(geo['vza'][iscan,ipxl]*np.pi/180.0) 
            results['amfGeo'][iscan,ipxl] = amfGeo
            results['NO2Geo'][iscan,ipxl] = results['NO2'][iscan,ipxl] / amfGeo
            results['NO2GeoSigma'][iscan,ipxl] = results['NO2Sigma'][iscan,ipxl] / amfGeo

        
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


            if cfg['export_spectra']:

                results['wvl'][iscan,ipxl,:len(wvl)] = wvl
                results['R_meas'][iscan,ipxl,:len(wvl)] = RMeas*szaFactor
                results['R_model'][iscan,ipxl,:len(wvl)] = RModel*szaFactor
                results['R_res'][iscan,ipxl,:len(wvl)] = RRes*szaFactor
                    

            # print(time.time()-start)
            results['errorFlag'][iscan,ipxl] = 0

        logger.info('Processed scanline {}/{} in {}s'.format(iscan-scanBeg+1,scanEnd-scanBeg+1,np.round((time.time() - startTimeScanline),2) ))

        return
    # end of scanline loop

    # multi-threading using ThreadPool
    if cfg['threads'] > 1:

        pool = multiprocessing.pool.ThreadPool(cfg['threads'])
        logger.info(f"Processing with {cfg['threads']} threads")

        # handle errors inside pool
        def handle_error(error):
            logger.error(error)

        # loop over scanlines in function for parallelization    
        for iscan in range(scanBeg,scanEnd+1,1):
            pool.apply_async(process_scanline,
                    args=(iscan,results),
                    error_callback=handle_error)

        pool.close()
        pool.join()
    
    # single-threaded
    else:
        for iscan in range(scanBeg,scanEnd+1,1):
            process_scanline(iscan,results)


    # H)  Finishing up

    # write results to output file
    libDOAS.writeOutput(cfg['io']['l2'],cfg,parameterNames,results,geo)


    logger.info(f"Output witten to {cfg['io']['l2']}")


    logger.info(f'IFDOE calculation finished in {np.round(time.time()-startTime,1)} s')
    

    return results


def amf_run(logger,cfg):

    startTime = time.time()

    logger.info(f"Reading DOAS results from L2 file: {cfg['io']['l2']}")
    doas = libAMF.read_doas(cfg['io']['l2'])

    logger.info(f"Reading atm file: {cfg['io']['sgm_atm']}")
    atm = libAMF.read_atm(cfg['io']['sgm_atm'])

    logger.info('Calculating AMF')
    amf_results = libAMF.get_amf(cfg['amf'], doas, atm)


    logger.info(f"Writing AMF results to: {cfg['io']['l2']}")
    libAMF.write_amf(cfg,amf_results)

    logger.info(f'AMF calculation finished in {np.round(time.time()-startTime,1)} s')

    return amf_results


def l1bl2_no2(logger,cfg):

    startTime = time.time()

    # use irradiance file from SGM. optional convolving
    if cfg['doas']['irr_from_sgm']:
        if cfg['doas']['convolve_irr']:
            convolved_irr_file = conv_irr(cfg['io']['sgm_rad'],cfg['isrf']['fwhm'])
            cfg['doas']['irr_file'] = convolved_irr_file
        else:
            cfg['doas']['irr_file'] = cfg['io']['sgm_rad']
    
    # optionally use radiance from SGM
    if cfg['doas']['rad_from_sgm']:
        cfg['doas']['rad_file'] = cfg['io']['sgm_rad']
    else:
        cfg['doas']['rad_file'] = cfg['io']['l1b']

    cfg['doas']['l2_file'] = cfg['io']['l2']
    cfg['doas']['gm_file'] = cfg['io']['gm']            


    # Python parallises internally with numpy, for single thread optimum is 4 numpy threads
    # for multi-threading use only 1 numpy thread, otherwise slow-down

    if cfg['doas']['threads'] == 1:
        numpy_cpu = 4
    else:
        numpy_cpu = 1

    with threadpool_limits(limits=numpy_cpu, user_api='blas'):
        doas_results = ifdoe_run(logger, cfg['doas'])
                        
    amf_results = amf_run(logger, cfg)

    logger.info(f'L1L2 calculation finished in {np.round(time.time()-startTime,1)} s')

    return

if __name__ == '__main__':
    

    # call with:
    # python l1bl2_no2.py l1bl2_no2.yaml

    # or with logging to file:
    # python l1bl2_no2.py l1bl2_no2.yaml l1bl2_no2.log

    # reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))

    if cfg['doas']['debug']['log']:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    # setup the logging to screen and to file
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

    if len(sys.argv) > 2:
        fh = logging.FileHandler(sys.argv[2], mode='w')
        fh.setLevel(logging.ERROR)
        fh.setFormatter(formatter)

        ch = logging.StreamHandler()
        ch.setLevel(loglevel)
        ch.setFormatter(formatter) 

        logging.basicConfig(level=loglevel, handlers = [ch,fh])

        logging.info(f'Logging to file: {sys.argv[2]}')

    else:
        ch = logging.StreamHandler()
        ch.setLevel(loglevel)
        ch.setFormatter(formatter) 
        logging.basicConfig(level=loglevel,handlers = [ch])
    
    logger = logging.getLogger()

    # Python parallises internally with numpy, for single thread optimum is 4 numpy threads
    # for multi-threading use only 1 numpy thread, otherwise slow-down

    if cfg['doas']['threads'] == 1:
        numpy_cpu = 4
    else:
        numpy_cpu = 1

    with threadpool_limits(limits=numpy_cpu, user_api='blas'):

        if cfg['doas']['run']:
            doas_results = ifdoe_run(logger, cfg)
                        
        if cfg['amf']['run']:
            amf_results = amf_run(logger,cfg)

