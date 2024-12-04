import numpy as np
import os
import sys
import yaml
from threadpoolctl import threadpool_limits
import time

from teds import log
from teds.lib import libDOAS, libAMF, libCloud, libConv

def get_slice(cfg):
    if 'alt' in cfg:
        slice_alt = slice(cfg['alt']['start'],cfg['alt']['stop']+1)
    else:
        slice_alt = slice(0,None)
    if 'act' in cfg:
        slice_act = slice(cfg['act']['start'],cfg['act']['stop']+1)
    else:
        slice_act = slice(0,None)

    return slice_alt, slice_act


def cloud_run(cfg):

    startTime = time.time()
    
    slice_alt, slice_act = get_slice(cfg)

    log.info(f"Reading DOAS results from L2 file: {cfg['io']['l2']}")
    doas = libCloud.read_doas(cfg['io']['l2'], slice_alt, slice_act)

    log.info(f"Reading atm file: {cfg['io']['sgm_atm']}")
    atm = libCloud.read_atm(cfg['io']['sgm_atm'], slice_alt, slice_act)

    log.info('Calculating cloud fraction')
    cloud_results = libCloud.get_cloud_fraction(cfg, doas, atm)

    log.info('Calculating cloud parameters')
    cloud_results = libCloud.get_cloud_parameters(cfg, doas, atm, cloud_results)

    log.info(f"Writing cloud results to: {cfg['io']['l2']}")
    libCloud.write_cloud(cfg, cloud_results, slice_alt, slice_act)

    log.info(f'Cloud calculation finished in {np.round(time.time()-startTime,1)} s')
    return cloud_results

def amf_run(cfg):

    startTime = time.time()

    slice_alt, slice_act = get_slice(cfg)

    log.info(f"Reading DOAS results from L2 file: {cfg['io']['l2']}")
    doas = libAMF.read_doas(cfg['io']['l2'], slice_alt, slice_act)

    log.info(f"Reading atm file: {cfg['io']['sgm_atm']}")
    atm = libAMF.read_atm(cfg['io']['sgm_atm'], slice_alt, slice_act)

    if doas['lat'].shape != atm['latitude'].shape:
        log.error(f'L1B radiance shape {doas['lat'].shape} does not match SGM ATM shape {atm['latitude'].shape}. Specify ALT and/or ACT in l2 config')
        sys.exit()

    log.info(f"Reading cloud results from L2 file: {cfg['io']['l2']}")
    cloud = libAMF.read_cloud(cfg['io']['l2'], slice_alt, slice_act)

    log.info('Calculating AMF')
    amf_results = libAMF.get_amf(cfg, doas, atm, cloud)

    log.info(f"Writing AMF results to: {cfg['io']['l2']}")
    libAMF.write_amf(cfg, amf_results, slice_alt, slice_act)

    log.info(f'AMF calculation finished in {np.round(time.time()-startTime,1)} s')
    return amf_results


def l1bl2_no2(cfg):

    startTime = time.time()

    if cfg['run_doas']:

        if os.path.isfile(cfg['io']['l2']):
            log.warning(f'File {cfg['io']['l2']} already exists, removing')
            os.remove(cfg['io']['l2'])

        # use irradiance file from SGM. optional convolving with Gaussian or ISRF
        if cfg['irr_from_sgm']:
            if cfg['convolve_irr'] == 'Gaussian':
                cfg['io']['sgm_irr'] = libConv.conv_irr(cfg['io']['sgm_rad'], mode='Gaussian', fwhm=cfg['gaussian_fwhm'])

            elif cfg['convolve_irr'] == 'ISRF':
                cfg['io']['sgm_irr'] = libConv.conv_irr(cfg['io']['sgm_rad'], mode='ISRF', isrf_file=cfg['isrf_file'])

            else:
                cfg['io']['sgm_irr'] = cfg['io']['sgm_rad']

        # use radiance file from SGM. optional convolving with Gaussian or ISRF
        if cfg['rad_from_sgm']:
            if cfg['convolve_rad'] == 'Gaussian':
                cfg['io']['sgm_rad'] = libConv.conv_rad(cfg['io']['sgm_rad'], mode='Gaussian', fwhm=cfg['gaussian_fwhm'])

            elif cfg['convolve_rad'] == 'ISRF':
                cfg['io']['sgm_rad'] = libConv.conv_rad(cfg['io']['sgm_rad'], mode='ISRF', isrf_file=cfg['isrf_file'])

        # Python parallises internally with numpy, for single thread optimum is 4 numpy threads
        # for multi-threading use only 1 numpy thread, otherwise slow-down

        if cfg['threads'] == 1:
            numpy_cpu = 4
        else:
            numpy_cpu = 1

        with threadpool_limits(limits=numpy_cpu, user_api='blas'):

            if cfg['retrieve']['no2']:
                doas_results_no2 = libDOAS.ifdoe_run(cfg, mode='no2')
            
            if cfg['retrieve']['o2o2']:
                doas_results_o2o2 = libDOAS.ifdoe_run(cfg, mode='o2o2')
    
    if cfg['run_clouds']:
        cloud_results = cloud_run(cfg)

    if cfg['run_amf']:               
        amf_results = amf_run(cfg)

    log.info(f'L1L2 calculation finished in {np.round(time.time()-startTime,1)} s')

    return

if __name__ == '__main__':
    

    # call with:
    # python l1bl2_no2.py l1bl2_no2.yaml

    # reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))


    # Python parallises internally with numpy, for single thread optimum is 4 numpy threads
    # for multi-threading use only 1 numpy thread, otherwise slow-down

    if cfg['doas']['threads'] == 1:
        numpy_cpu = 4
    else:
        numpy_cpu = 1

    with threadpool_limits(limits=numpy_cpu, user_api='blas'):

        if cfg['doas']['run']:
            doas_results = libDOAS.ifdoe_run(cfg)
                        
        if cfg['amf']['run']:
            amf_results = amf_run(cfg)

