// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_proc.h"
#include "ckd.h"
#include "l1a_file_metadata.h"
#include "binningtable.h"
#include "l1x.h"
#include "l1a.h"

// Constructor.
L1X::L1X( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
L1X::~L1X() {} // Destructor.

// Initialization of L1X file.
int L1X::init( // {{{
    l1x_t il1x, // L1X level identifier.
    size_t nframe, // Number of frames (assessed by L1A manager).
    string &filename, // L1X file name.
    L1A_file_metadata *file_meta, // L1A file metadata.
    CKD *ckd, // Calibration key data (for detector and spectrum dimensions).
    size_t largest_spectrum_size // NetCDF dimension size for extracted spectra (only for FOV and later L1X steps).
)
{

    step = il1x; // Save step.

    // Verify if CKD is mature enough to get to the desired level. This is important, because
    // we will get dimensions from the CKD, which may be uninitialized or zero and cause
    // all kinds of issues. Moreover, it is not necessary to define L1X of a higher level
    // than you can calculate. It would generate an empty file, which is not particularly
    // interesting.
    if (step == L1X_RAW) {
        // Probably, this assertion will never take place. If there is any L1A, the CKD must
        // have passed the 'dim' step. If there is no L1A, there is no L1X. You can select
        // L1X by defining the equal number of files, which is zero. That zero will be
        // recognized as 'no L1X'.
        check_error(ckd->lev <= LEVEL_DIMCAL,"Error: CKD must have passed the 'dim' level to be able to do anything, including writing raw images. This has to do with the binning table and reading routine.");
    }
    if (step == L1X_DARK) {
        check_error(ckd->lev <= LEVEL_DARKCAL,"Error: CKD must have passed the 'dark' level to be able to write L1X after dark correction.");
    }
    if (step == L1X_NOISE) {
        check_error(ckd->lev <= LEVEL_NOISECAL,"Error: CKD must have passed the 'noise' level to be able to write L1X after noise estimation.");
    }
    if (step == L1X_NONLIN) {
        check_error(ckd->lev <= LEVEL_NONLINCAL,"Error: CKD must have passed the 'nonlin' level to be able to write L1X after non-linearity correction.");
    }
    if (step == L1X_PRNU) {
        check_error(ckd->lev <= LEVEL_PRNUCAL,"Error: CKD must have passed the 'prnu' level to be able to write L1X after PRNU correction.");
    }
    if (step == L1X_UNBIN || step == L1X_STRAY || step == L1X_REBIN) {
        check_error(ckd->lev <= LEVEL_STRAYCAL,"Error: CKD must have passed the 'stray' level to be able to write L1X after unbinning, straylight correction or rebinning.");
    }
    if (step == L1X_FOV) {
        check_error(ckd->lev <= LEVEL_FOVCAL,"Error: CKD must have passed the 'fov' level to be able to write L1X after spectrum extraction.");
    }
    if (step == L1X_RAD) {
        check_error(ckd->lev <= LEVEL_RADCAL,"Error: CKD must have passed the 'rad' level to be able to write L1X after radiometric calibration.");
    }

    // Set execution flag.
    execute = true;

    // Create the L1X output file.
    nc = make_unique<NetCDF_object>(this);
    handle(nc->open(filename,NcFile::replace));

    // Dimensions.
    NcDim dimid_nframe;
    netcdf_check(nc,dimid_nframe = nc->ncid->addDim("number_of_images",nframe));
    NcDim dimid_npix;
    netcdf_check(nc,dimid_npix = nc->ncid->addDim("samples_per_image",file_meta->npix));
    // Interpolated contents is already interpolated so it has the same
    // dimension as the frames. TODO: Skip interpolation when reading?
    NcDim dimid_nframe_eng;
    netcdf_check(nc,dimid_nframe_eng = nc->ncid->addDim("hk_packets",nframe));
    // Dimension SC_hkt_block is not used.
    NcDim dimid_nframe_nav;
    netcdf_check(nc,dimid_nframe_nav = nc->ncid->addDim("SC_records",file_meta->l1a_navigation?nframe:0));
    // Trivial geometric dimensions, used in navigation data.
    NcDim dimid_vec;
    netcdf_check(nc,dimid_vec = nc->ncid->addDim("vector_elements",DIM_VEC));
    NcDim dimid_quat;
    netcdf_check(nc,dimid_quat = nc->ncid->addDim("quaternion_elements",DIM_QUAT));

    // Write global attributes that are read in L1A file metadata.
    if (file_meta->startDirection.compare("Unknown") != 0) {
        netcdf_check(nc,nc->ncid->putAtt("startDirection",file_meta->startDirection));
    }
    if (file_meta->endDirection.compare("Unknown") != 0) {
        netcdf_check(nc,nc->ncid->putAtt("endDirection",file_meta->endDirection));
    }
    if (file_meta->orbit_number != NC_FILL_INT64) {
        netcdf_check(nc,nc->ncid->putAtt("orbit_number",ncInt64,file_meta->orbit_number));
    }
    if (file_meta->time_coverage_start.compare("Unknown") != 0) {
        netcdf_check(nc,nc->ncid->putAtt("time_coverage_start",file_meta->time_coverage_start));
    }
    if (file_meta->time_coverage_stop.compare("Unknown") != 0) {
        netcdf_check(nc,nc->ncid->putAtt("time_coverage_stop",file_meta->time_coverage_stop));
    }

    // Write L1X maturity signature in global attributes.
    netcdf_check(nc,nc->ncid->putAtt("l1x_maturity_level",ncInt,step));

    // Image attributes.
    NcGroup grp; // Only used when creating variables. After that, the object can be overwritten.
    netcdf_check(nc,grp = nc->ncid->addGroup("image_attributes"));
    netcdf_check(nc,var_image_time = grp.addVar("image_time",ncDouble,dimid_nframe));
    netcdf_check(nc,var_image_time.putAtt("long_name","image time (seconds of day)"));
    netcdf_check(nc,var_image_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc,var_image_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc,var_image_time.putAtt("units","seconds"));
    if (file_meta->time_reference.compare("Unknown") != 0) {
        netcdf_check(nc,var_image_time.putAtt("reference",file_meta->time_reference));
    }
    netcdf_check(nc,var_image_CCSDS_sec = grp.addVar("image_CCSDS_sec",ncUint,dimid_nframe));
    netcdf_check(nc,var_image_CCSDS_sec.putAtt("long_name","image CCSDS time (seconds since 1970)"));
    netcdf_check(nc,var_image_CCSDS_sec.putAtt("valid_min",ncUint,1900000000));
    netcdf_check(nc,var_image_CCSDS_sec.putAtt("valid_max",ncUint,2400000000));
    netcdf_check(nc,var_image_CCSDS_sec.putAtt("units","seconds"));
    netcdf_check(nc,var_image_CCSDS_usec = grp.addVar("image_CCSDS_usec",ncInt,dimid_nframe));
    netcdf_check(nc,var_image_CCSDS_usec.putAtt("long_name","image CCSDS time (microseconds)"));
    netcdf_check(nc,var_image_CCSDS_usec.putAtt("valid_min",ncInt,0));
    netcdf_check(nc,var_image_CCSDS_usec.putAtt("valid_max",ncInt,999999));
    netcdf_check(nc,var_image_CCSDS_usec.putAtt("units","microseconds"));
    netcdf_check(nc,var_binning_table = grp.addVar("binning_table",ncUbyte,dimid_nframe));
    netcdf_check(nc,var_binning_table.putAtt("long_name","binning-table ID"));
    netcdf_check(nc,var_binning_table.putAtt("valid_min",ncUbyte,0));
    netcdf_check(nc,var_binning_table.putAtt("valid_max",ncUbyte,255));
    netcdf_check(nc,var_digital_offset = grp.addVar("digital_offset",ncShort,dimid_nframe));
    netcdf_check(nc,var_digital_offset.putAtt("long_name","digital offset"));
    netcdf_check(nc,var_digital_offset.putAtt("units","1"));
    netcdf_check(nc,var_nr_coadditions = grp.addVar("nr_coadditions",ncUshort,dimid_nframe));
    netcdf_check(nc,var_nr_coadditions.putAtt("long_name","number of coadditions"));
    netcdf_check(nc,var_nr_coadditions.putAtt("units","1"));
    netcdf_check(nc,var_exposure_time = grp.addVar("exposure_time",ncDouble,dimid_nframe));
    netcdf_check(nc,var_exposure_time.putAtt("long_name","exposure time"));
    netcdf_check(nc,var_exposure_time.putAtt("units","s"));

    // Engineering data.
    netcdf_check(nc,grp = nc->ncid->addGroup("engineering_data"));
    netcdf_check(nc,var_hk_tlm_time = grp.addVar("HK_tlm_time",ncDouble,dimid_nframe_eng));
    netcdf_check(nc,var_hk_tlm_time.putAtt("long_name","HK telemetry packet time (seconds of day)"));
    netcdf_check(nc,var_hk_tlm_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc,var_hk_tlm_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc,var_hk_tlm_time.putAtt("units","seconds"));
    netcdf_check(nc,var_temp_detector = grp.addVar("temp_detector",ncDouble,dimid_nframe_eng));
    netcdf_check(nc,var_temp_detector.putAtt("long_name","Detector temperature"));
    netcdf_check(nc,var_temp_detector.putAtt("valid_min",ncDouble,260.0));
    netcdf_check(nc,var_temp_detector.putAtt("valid_max",ncDouble,300.0));
    netcdf_check(nc,var_temp_detector.putAtt("units","K"));

    // Science data. This is where the shape differs per L1X level.
    netcdf_check(nc,grp = nc->ncid->addGroup("science_data"));
    // TODO: Do this. The reason we do this in the middle is that this
    // is the group order in the L1A document.

    // For any time we define a mask and want to write the attribute.
    const vector<uint8_t> mask_values = {0,1};

    // The (binned) image is there during detector calibration,
    // temporarily not during straylight.
    if (
        step == L1X_RAW ||
        step == L1X_DARK ||
        step == L1X_NOISE ||
        step == L1X_NONLIN ||
        step == L1X_PRNU ||
        step == L1X_REBIN
    ) {
        // Raw image.
        netcdf_check(nc,var_img = grp.addVar("detector_images",ncDouble,{dimid_nframe,dimid_npix}));
        netcdf_check(nc,var_img.putAtt("long_name","Detector image at current stage."));
        netcdf_check(nc,var_img.putAtt("units","counts"));
        // Pixel mask.
        netcdf_check(nc,var_img_mask = grp.addVar("detector_images_mask",ncUbyte,{dimid_nframe,dimid_npix}));
        netcdf_check(nc,var_img.putAtt("long_name","Detector pixel mask."));
        netcdf_check(nc,var_img_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
        netcdf_check(nc,var_img_mask.putAtt("flag_meanings","Good Bad"));
    }

    // The noise starts existing from noise. Furthermore, it stays in
    // the binned perspective during unbinning and straylight.
    if (
        step == L1X_NOISE ||
        step == L1X_NONLIN ||
        step == L1X_PRNU ||
        step == L1X_UNBIN ||
        step == L1X_STRAY ||
        step == L1X_REBIN
    ) {
        // Noise estimate.
        netcdf_check(nc,var_img_noise = grp.addVar("detector_images_noise",ncDouble,{dimid_nframe,dimid_npix}));
        netcdf_check(nc,var_img_noise.putAtt("long_name","Noise estimate of image corrected for dark offset and current."));
        netcdf_check(nc,var_img_noise.putAtt("units","counts"));
    }

    // Special case, dark offset corrected, but not yet the current.
    if (step == L1X_DARK) {
        // Image corrected for just offset.
        netcdf_check(nc,var_img_with_current = grp.addVar("detector_images_with_current",ncDouble,{dimid_nframe,dimid_npix}));
        netcdf_check(nc,var_img_with_current.putAtt("long_name","Image before correcting dark current."));
        netcdf_check(nc,var_img_with_current.putAtt("units","counts"));
    }

    // Unbinned images.
    if (
        step == L1X_UNBIN ||
        step == L1X_STRAY
    ) {
        // 2D detector image, spatial and spectral detector dimension needed as well.
        NcDim dimid_detector_spec;
        netcdf_check(nc,dimid_detector_spec = grp.addDim("detector_spec",ckd->dim_detector_spec));
        NcDim dimid_detector_spat;
        netcdf_check(nc,dimid_detector_spat = grp.addDim("detector_spat",ckd->dim_detector_spat));
        // 2D detector image. Still we use the same variable idenitifer, but we will
        // name this differently so that it is distinguishable from the binned image.
        netcdf_check(nc,var_img = grp.addVar("detector_images_2d",ncDouble,{dimid_nframe,dimid_detector_spat,dimid_detector_spec}));
        netcdf_check(nc,var_img.putAtt("long_name","Two-dimensional detector image at current stage."));
        netcdf_check(nc,var_img.putAtt("units","counts"));
        // No noise. There is no such thing as 2d noise. It would be
        // correlated rubbish.
        // Mask exists to aid the straylight correction. TODO: In the new binning protocol,
        // the unbinned mask will be thrown away after straylight and the binned one is
        // dug up again. So, for L1X, the binned noise and mask should survive until the
        // end of the straylight correction.
        netcdf_check(nc,var_img_mask = grp.addVar("detector_images_2d_mask",ncUbyte,{dimid_nframe,dimid_detector_spat,dimid_detector_spec}));
        netcdf_check(nc,var_img_mask.putAtt("long_name","Two-dimensional detector pixel mask."));
        netcdf_check(nc,var_img_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
        netcdf_check(nc,var_img_mask.putAtt("flag_meanings","Good Bad"));
    }

    // Spectra.
    if (
        step == L1X_FOV ||
        step == L1X_RAD
    ) {
        // Extracted spectra, noise and mask. The spectral size is the
        // largest size that exists over the entire L1A ensemble, not just
        // the L1A that are in this file. In practical cases, this is never
        // a sacrifice.
        NcDim dimid_spectra_max;
        NcDim dimid_fov;
        NcDim dimid_pol;
        netcdf_check(nc,dimid_spectra_max = grp.addDim("spectral_samples_per_spectrum",largest_spectrum_size));
        netcdf_check(nc,dimid_fov = grp.addDim("spatial_samples_per_image",ckd->dim_fov));
        netcdf_check(nc,dimid_pol = grp.addDim("polarization_state",DIM_POL));

        // Actual size of the spectra. This depends on FOV (actually viewport),
        // but also on which L1A you are, because different L1A can have different
        // binning tables.
        netcdf_check(nc,var_spectra_sizes = grp.addVar("spectra_sizes",ncUint,{dimid_nframe,dimid_fov}));
        netcdf_check(nc,var_spectra_sizes.putAtt("long_name","Number of spectral bins per extracted spectrum"));

        // Unit of the spectrum depends on the step.
        string unit = step == L1X_FOV?"counts":"W.m-2.sr-1.um-1";

        // Spectra.
        netcdf_check(nc,var_spectra = grp.addVar("spectra",ncDouble,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
        netcdf_check(nc,var_spectra.putAtt("long_name","Extracted spectra at current stage."));
        netcdf_check(nc,var_spectra.putAtt("units",unit));
        // Noise.
        netcdf_check(nc,var_spectra_noise = grp.addVar("spectra_noise",ncDouble,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
        netcdf_check(nc,var_spectra_noise.putAtt("long_name","Noise estimate on extracted spectra at current stage."));
        netcdf_check(nc,var_spectra_noise.putAtt("units",unit));
        // Mask.
        netcdf_check(nc,var_spectra_mask = grp.addVar("spectra_mask",ncUbyte,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
        netcdf_check(nc,var_spectra_mask.putAtt("long_name","Pixel mask on extracted spectra."));
        netcdf_check(nc,var_spectra_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
        netcdf_check(nc,var_spectra_mask.putAtt("flag_meanings","Good Bad"));

        // Additional spectrum for RAD, the interpolated one.
        if (step == L1X_RAD) {
            // Repeat joke after wavelength interpolation.
            // Spectra.
            netcdf_check(nc,var_spectra_target = grp.addVar("spectra_target",ncDouble,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
            netcdf_check(nc,var_spectra_target.putAtt("long_name","Extracted spectra at current stage after wavelength interpolation."));
            netcdf_check(nc,var_spectra_target.putAtt("units",unit));
            // Noise.
            netcdf_check(nc,var_spectra_target_noise = grp.addVar("spectra_target_noise",ncDouble,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
            netcdf_check(nc,var_spectra_target_noise.putAtt("long_name","Noise estimate on extracted spectra at current stage after wavelength interpolation."));
            netcdf_check(nc,var_spectra_target_noise.putAtt("units",unit));
            // Mask.
            netcdf_check(nc,var_spectra_target_mask = grp.addVar("spectra_target_mask",ncUbyte,{dimid_nframe,dimid_fov,dimid_pol,dimid_spectra_max}));
            netcdf_check(nc,var_spectra_target_mask.putAtt("long_name","Pixel mask on extracted spectra after wavelength interpolation."));
            netcdf_check(nc,var_spectra_target_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc,var_spectra_target_mask.putAtt("flag_meanings","Good Bad"));
        }
    }

    // Navigation data.
    netcdf_check(nc,grp = nc->ncid->addGroup("navigation_data"));
    netcdf_check(nc,var_att_time = grp.addVar("att_time",ncDouble,dimid_nframe_nav));
    netcdf_check(nc,var_att_time.putAtt("long_name","Attitude sample time (seconds of day)"));
    netcdf_check(nc,var_att_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc,var_att_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc,var_att_time.putAtt("units","seconds"));
    netcdf_check(nc,var_att_qaut = grp.addVar("att_quat",ncDouble,{dimid_nframe_nav,dimid_quat}));
    netcdf_check(nc,var_att_qaut.putAtt("long_name","Attitude quaternions (J2000 to spacecraft)"));
    netcdf_check(nc,var_att_qaut.putAtt("valid_min",ncDouble,-1.0));
    netcdf_check(nc,var_att_qaut.putAtt("valid_max",ncDouble,1.0));
    netcdf_check(nc,var_att_qaut.putAtt("units","1"));
    netcdf_check(nc,var_orb_time = grp.addVar("orb_time",ncDouble,dimid_nframe_nav));
    netcdf_check(nc,var_orb_time.putAtt("long_name","Orbit vector time (seconds of day)"));
    netcdf_check(nc,var_orb_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc,var_orb_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc,var_orb_time.putAtt("units","seconds"));
    netcdf_check(nc,var_orb_pos = grp.addVar("orb_pos",ncDouble,{dimid_nframe_nav,dimid_vec}));
    netcdf_check(nc,var_orb_pos.putAtt("long_name","Orbit positions vectors (ECR)"));
    netcdf_check(nc,var_orb_pos.putAtt("valid_min",ncDouble,-7200000.0));
    netcdf_check(nc,var_orb_pos.putAtt("valid_max",ncDouble,7200000.0));
    netcdf_check(nc,var_orb_pos.putAtt("units","meters"));
    netcdf_check(nc,var_orb_vel = grp.addVar("orb_vel",ncDouble,{dimid_nframe_nav,dimid_vec}));
    netcdf_check(nc,var_orb_vel.putAtt("long_name","Orbit velocity vectors (ECR)"));
    netcdf_check(nc,var_orb_vel.putAtt("valid_min",ncDouble,-7600.0));
    netcdf_check(nc,var_orb_vel.putAtt("valid_max",ncDouble,7600.0));
    netcdf_check(nc,var_orb_vel.putAtt("units","meters/seconds"));

    // GSE-data, all in file_meta.
    // We will not write fillvalue data, though we make the group even
    // if it ends up being empty.
    netcdf_check(nc,grp = nc->ncid->addGroup("gse_data"));
    NcVar var;
    if (file_meta->viewport != 0) {
        netcdf_check(nc,var = grp.addVar("viewport",ncUbyte));
        netcdf_check(nc,var.putAtt("long_name","viewport status"));
        vector<uint8_t> rng = {0,16};
        netcdf_check(nc,var.putAtt("valid_range",ncUbyte,2,rng.data()));
        netcdf_check(nc,var.putAtt("comment","bitmask: 1, 2, 4, 8, 16"));
        netcdf_check(nc,var.putVar(&file_meta->viewport));
    }
    if (file_meta->dim_refspec != 0) {
        NcDim dimid_refspec;
        netcdf_check(nc,dimid_refspec = grp.addDim("wavelength",file_meta->dim_refspec));
        // If the dimension exists, we will write the reference spectrum.
        // On the input side, it was legal to define the dimension but not
        // the spectrum in which case fillvalues are written. It is not
        // worth the effort checking if the entire spectrum is fillvalue
        // and then not write the spectrum in the highly unlikely event that
        // someone without a spectrum actually writes down a nonzero dimension.
        NcVar var; // Variable to be written directly.
        netcdf_check(nc,var = grp.addVar("wavelength",ncDouble,dimid_refspec));
        netcdf_check(nc,var.putAtt("long_name","wavelength of stimulus"));
        netcdf_check(nc,var.putAtt("units","nm"));
        netcdf_check(nc,var.putVar(file_meta->refspec_wavelength.data()));
        netcdf_check(nc,var = grp.addVar("signal",ncDouble,dimid_refspec));
        netcdf_check(nc,var.putAtt("long_name","signal of stimulus"));
        netcdf_check(nc,var.putAtt("units","W/(nm.m^2)"));
        netcdf_check(nc,var.putVar(file_meta->refspec_radiance.data()));
    }
    if (file_meta->illumination_level != NC_FILL_DOUBLE) {
        netcdf_check(nc,grp.putAtt("Illumination_level",ncDouble,file_meta->illumination_level));
    }
    if (file_meta->refspec_dolp != NC_FILL_DOUBLE) {
        netcdf_check(nc,grp.putAtt("DoLP",ncDouble,file_meta->refspec_dolp));
    }
    if (file_meta->refspec_aolp != NC_FILL_DOUBLE) {
        netcdf_check(nc,grp.putAtt("AoLP",ncDouble,file_meta->refspec_aolp/DEGREES));
    }
    if (file_meta->act_angle != NC_FILL_DOUBLE) {
        netcdf_check(nc,grp.putAtt("ACT_rotationAngle",ncDouble,file_meta->act_angle/DEGREES));
    }
    if (file_meta->alt_angle != NC_FILL_DOUBLE) {
        netcdf_check(nc,grp.putAtt("ALT_rotationAngle",ncDouble,file_meta->alt_angle/DEGREES));
    }

    return 0;

} // }}}

int L1X::write_metadata( // {{{
    L1A *l1a
)
{

    // Only write if you exist.
    if (execute) {

        // This can be done independent of L1X level. All variables filled here are in groups
        // different from science_data.
        const vector<size_t> strt = {l1a->l1x_iframe};
        const vector<size_t> cnt = {1};

        // Image attributes.
        netcdf_check(nc,var_image_time.putVar(strt,cnt,&l1a->time));
        // Warning for inconsistency. The absolute time is averaged between first and last
        // co-added frame, which can be a half microsecond. This half microsecond will be
        // rounded down now.
        netcdf_check(nc,var_image_CCSDS_sec.putVar(strt,cnt,&l1a->seconds));
        // Convert secondfraction to microseconds.
        int microseconds = static_cast<int>(1.0e6 * l1a->secondfraction + 0.5); // +0.5 to round to nearest integer, preventing numeric misery.
        netcdf_check(nc,var_image_CCSDS_usec.putVar(strt,cnt,&microseconds));
        netcdf_check(nc,var_binning_table.putVar(strt,cnt,&l1a->binning_table_id));
        netcdf_check(nc,var_nr_coadditions.putVar(strt,cnt,&l1a->nr_coadditions));
        netcdf_check(nc,var_exposure_time.putVar(strt,cnt,&l1a->exposure_time));

        // Engineering data.
        // Time axis is equal to the image time, because this is after interpolation.
        netcdf_check(nc,var_hk_tlm_time.putVar(strt,cnt,&l1a->time));

        // Science data is not metadata.

        // Navigation data.
        if (l1a->l1a_navigation) {
            // Here, also all interpolations have already been done.
            netcdf_check(nc,var_att_time.putVar(strt,cnt,&l1a->time));
            // All vectors and quaternions are extracted first.
            const vector<size_t> strt_quat = {l1a->l1x_iframe,0};
            const vector<size_t> cnt_quat = {1,DIM_QUAT};
            vector<double> extracted_quaternion(DIM_QUAT);
            l1a->attitude_quaternion.get(extracted_quaternion.data());
            netcdf_check(nc,var_att_qaut.putVar(strt_quat,cnt_quat,extracted_quaternion.data()));
            netcdf_check(nc,var_orb_time.putVar(strt,cnt,&l1a->time));
            // Recycling a four-element array as three-element array is ugly.
            vector<double> extracted_vector(DIM_VEC);
            l1a->satellite_position.get(extracted_vector.data());
            // Recycling another one's start vector is ugly.
            const vector<size_t> strt_vec = {l1a->l1x_iframe,0};
            const vector<size_t> cnt_vec = {1,DIM_VEC};
            netcdf_check(nc,var_orb_pos.putVar(strt_vec,cnt_vec,extracted_vector.data()));
            // Recycling a vector for another vector is okay as long as the name just says vector.
            l1a->satellite_velocity.get(extracted_vector.data());
            // Same for start and count.
            netcdf_check(nc,var_orb_vel.putVar(strt_vec,cnt_vec,extracted_vector.data()));
        }

        // GSE data is already filled.

    }

    return 0;

} // }}}

int L1X::write( // {{{
    L1A *l1a, // L1A instance that contains the information.
    size_t ifov // Field-of-view coordinate (only for spectra, L1X_FOV and L1X_RAD).
)
{

    if (execute) {
        // Write binned images, relevant for all these steps.
        if (
            step == L1X_RAW ||
            step == L1X_DARK ||
            step == L1X_NOISE ||
            step == L1X_NONLIN ||
            step == L1X_PRNU ||
            step == L1X_UNBIN ||
            step == L1X_STRAY ||
            step == L1X_REBIN
        ) {
            // Define start and count vectors, relevant for all the
            // binned images.
            size_t npix_frac = l1a->ipix_end-l1a->ipix_start;
            vector<size_t> strt = {l1a->l1x_iframe,l1a->ipix_start};
            vector<size_t> cnt = {1,npix_frac};

            // Image and mask, not for unbinning and straylight.
            if (
                step == L1X_RAW ||
                step == L1X_DARK ||
                step == L1X_NOISE ||
                step == L1X_NONLIN ||
                step == L1X_PRNU ||
                step == L1X_REBIN
            ) {

                // The image.
                netcdf_check(nc,var_img.putVar(strt,cnt,&l1a->image[l1a->ipix_start]));
                // Mask.
                // Even though the mask is officially defined for the entire
                // (binned) detector, we only write the fraction relevant for
                // this L1A.
                vector<uint8_t> writemaskvector(npix_frac);
                uint8_t *writemask = writemaskvector.data() - l1a->ipix_start;
                for (size_t ipix=l1a->ipix_start ; ipix<l1a->ipix_end ; ipix++) writemask[ipix] = l1a->pixelmask[ipix]?1:0;
                netcdf_check(nc,var_img_mask.putVar(strt,cnt,&writemask[l1a->ipix_start]));

            }
            // Image noise, only relevant from noise step onwards. And
            // also relevant for unbinning and straylight, because the
            // binned noise image remains in these steps.
            if (
                step == L1X_NOISE ||
                step == L1X_NONLIN ||
                step == L1X_PRNU ||
                step == L1X_UNBIN ||
                step == L1X_STRAY ||
                step == L1X_REBIN
            ) {
                netcdf_check(nc,var_img_noise.putVar(strt,cnt,&l1a->noise[l1a->ipix_start]));
            }

            // For dark, also write the image with dark current not yet subtracted.
            if (step == L1X_DARK) {
                netcdf_check(nc,var_img_with_current.putVar(strt,cnt,&l1a->image_with_dark_current[l1a->ipix_start]));
            }
        }

        // Write unbinned image, noise and mask.
        if (
            step == L1X_UNBIN ||
            step == L1X_STRAY
        ) {
            // Short access to detector dimensions.
            size_t &dim_detector_spec = l1a->bin->binned_ckd->dim_detector_spec;
            size_t &dim_detector_spat = l1a->bin->binned_ckd->dim_detector_spat;
            // Ironically, we have the binning table for the unbinned
            // detector dimensions. It contains CKD and 2D detector
            // dimensions are always unbinned.
            vector<size_t> strt = {l1a->l1x_iframe,0,0};
            vector<size_t> cnt = {1,dim_detector_spat,dim_detector_spec};
            // In the case there is a fractional image, we do not know
            // the fraction of this L1A. We do, but we are too lazy to
            // consult the entire binning table for this. Instead, we
            // will just look what it already written, and overwrite
            // everything that is not a fillvaule in the current L1A.
            vector<double> frame(l1a->bin->npix_unbinned); // Recycled for image and noise.
            // Image.
            netcdf_check(nc,var_img.getVar(strt,cnt,frame.data()));
            // Fill values that are defined by the frames to be written.
            for (size_t ipix=0 ; ipix<l1a->bin->npix_unbinned ; ipix++) {
                double &el = l1a->image[ipix];
                if (el != NC_FILL_DOUBLE) frame[ipix] = el;
            }
            // Write the field as it is now.
            netcdf_check(nc,var_img.putVar(strt,cnt,frame.data()));

            // For the mask, it works the same, but the additional feature
            // that in the case of a fraction, excluded pixels are flagged
            // as dead and should not overwrite living pixels from another
            // fraction.
            vector<uint8_t> agg_mask(l1a->bin->npix_unbinned);
            netcdf_check(nc,var_img_mask.getVar(strt,cnt,agg_mask.data()));
            // If any fraction says that the pixel lives, it lives.
            // Dead pixels are dead. We do not want fill values.
            // It is better to be dead than to be a fill value.
            for (size_t ipix=0 ; ipix<l1a->bin->npix_unbinned ; ipix++) {
                if (agg_mask[ipix] == NC_FILL_UBYTE) agg_mask[ipix] = l1a->pixelmask[ipix]?1:0; // Overrule fillvalue with living or dead.
                else if (!l1a->pixelmask[ipix]) agg_mask[ipix] = 0; // Overrule dead with living.
            }
            // Write the field as it is now.
            netcdf_check(nc,var_img_mask.putVar(strt,cnt,agg_mask.data()));
        }

        // Spectra. This requires the ifov argument.
        if (
            step == L1X_FOV ||
            step == L1X_RAD
        ) {
            size_t &dim_spec = l1a->bin->binned_ckd->fov_dims_spec[ifov];
            check_error(ifov == NC_FILL_UINT64,"Error: Cannot write L1X spectra without field-of-view index.");

            // The size of the spectra, one write instruction with different start and count vectors.
            netcdf_check(nc,var_spectra_sizes.putVar({l1a->l1x_iframe,ifov},{1,1},&dim_spec));

            // Start and count vectors for the serious vectors.
            vector<size_t> strt = {l1a->l1x_iframe,ifov,0,0};
            vector<size_t> cnt = {1,1,DIM_POL,dim_spec};

            // Spectra, the signal.
            netcdf_check(nc,var_spectra.putVar(strt,cnt,l1a->l1x_specs->signal.data()));
            // Noise.
            netcdf_check(nc,var_spectra_noise.putVar(strt,cnt,l1a->l1x_specs->noise.data()));
            // Mask. No issues about incomplete spectra, just the conversion from a boolean to an unsigned byte.
            vector<uint8_t> writemask(DIM_POL*dim_spec);
            for (size_t iel=0 ; iel<DIM_POL*dim_spec ; iel++) writemask[iel] = l1a->l1x_specs->mask[iel]?1:0;
            netcdf_check(nc,var_spectra_mask.putVar(strt,cnt,writemask.data()));

            // Repeat joke for the 'target' spectra for RAD.
            if (step == L1X_RAD) {

                // Spectra, the signal.
                netcdf_check(nc,var_spectra_target.putVar(strt,cnt,l1a->l1x_specs_target->signal.data()));
                // Noise.
                netcdf_check(nc,var_spectra_target_noise.putVar(strt,cnt,l1a->l1x_specs_target->noise.data()));
                // Mask. No issues about incomplete spectra, just the conversion from a boolean to an unsigned byte.

                // Recycle the writemask.
                for (size_t iel=0 ; iel<DIM_POL*dim_spec ; iel++) writemask[iel] = l1a->l1x_specs_target->mask[iel]?1:0;
                netcdf_check(nc,var_spectra_target_mask.putVar(strt,cnt,writemask.data()));
            }
        }
    }

    return 0;

} // }}}

