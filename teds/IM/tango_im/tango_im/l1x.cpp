#include "header.h"
#include "netcdf_object.h"
#include "settings_noise.h"
#include "ckd.h"
#include "l1x_inputfile.h"
#include "frame.h"
#include "l1x.h"
#include <random>

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
    size_t nframe, // Number of frames.
    string &filename, // L1X file name.
    size_t a_npix_binned, // Number of pixels in binned format. Equal to entire detector for non-L1A.
    int a_binning_table_id, // Binning table ID. Zero for non-L1A.
    CKD *ckd, // Calibration key data (for detector and spectrum dimensions).
    L1X_inputfile *l1x_inputfile // L1X input file structure.
)
{

    step = il1x; // Save step.

    // Set execution flag.
    execute = true;

    // Save binning table ID.
    binning_table_id = a_binning_table_id;
    // Save binned number of pixels.
    npix_binned = a_npix_binned;
    npix = ckd->npix;

    // Create the L1X output file.
    nc = make_unique<NetCDF_object>(this);
    handle(nc->open(filename,NcFile::replace));

    // Dimensions.
    NcDim dimid_nframe;
    netcdf_check(nc,dimid_nframe = nc->ncid->addDim("number_of_images",nframe));
    NcDim dimid_npix;
    netcdf_check(nc,dimid_npix = nc->ncid->addDim("samples_per_image",npix_binned));

    // Write L1X maturity signature in global attributes. Not for L1A,
    // because it is not in the official L1A data format.
    if (step != L1X_L1A) {
        netcdf_check(nc,nc->ncid->putAtt("l1x_maturity_level",ncInt,step));
    }

    // Image attributes.
    NcGroup grp; // Only used when creating variables. After that, the object can be overwritten.
    netcdf_check(nc,grp = nc->ncid->addGroup("image_attributes"));
    netcdf_check(nc,var_image_time = grp.addVar("image_time",ncDouble,dimid_nframe));
    netcdf_check(nc,var_image_time.putAtt("long_name","image time (seconds of day)"));
    netcdf_check(nc,var_image_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc,var_image_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc,var_image_time.putAtt("units","seconds"));
    netcdf_check(nc,var_binning_table = grp.addVar("binning_table",ncUbyte,dimid_nframe));
    netcdf_check(nc,var_binning_table.putAtt("long_name","binning-table ID"));
    netcdf_check(nc,var_binning_table.putAtt("valid_min",ncUbyte,0));
    netcdf_check(nc,var_binning_table.putAtt("valid_max",ncUbyte,255));
    netcdf_check(nc,var_nr_coadditions = grp.addVar("nr_coadditions",ncUshort,dimid_nframe));
    netcdf_check(nc,var_nr_coadditions.putAtt("long_name","number of coadditions"));
    netcdf_check(nc,var_nr_coadditions.putAtt("units","1"));
    netcdf_check(nc,var_exposure_time = grp.addVar("exposure_time",ncDouble,dimid_nframe));
    netcdf_check(nc,var_exposure_time.putAtt("long_name","exposure time"));
    netcdf_check(nc,var_exposure_time.putAtt("units","s"));

    // Science data. This is where the shape differs per L1X level.
    netcdf_check(nc,grp = nc->ncid->addGroup("science_data"));

    // L1A output. This has a different data type for the image.
    if (step == L1X_L1A) {
        // Raw image.
        netcdf_check(nc,var_img = grp.addVar("detector_images",ncInt,{dimid_nframe,dimid_npix}));
        netcdf_check(nc,var_img.putAtt("long_name","Detector image"));
        netcdf_check(nc,var_img.putAtt("units","counts"));
    }

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
        NcDim dimid_spectra_max;
        NcDim dimid_fov;
        NcDim dimid_pol;
        netcdf_check(nc,dimid_spectra_max = grp.addDim("spectral_samples_per_spectrum",ckd->dim_detector_spec));
        netcdf_check(nc,dimid_fov = grp.addDim("spatial_samples_per_image",ckd->dim_fov));
        netcdf_check(nc,dimid_pol = grp.addDim("polarization_state",DIM_POL));

        // Actual size of the spectra.
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
    // No L1X_TRUTH, because that can never be output.

    return 0;

} // }}}

int L1X::write_metadata( // {{{
    Frame *frm
)
{

    // Only write if you exist.
    if (execute) {

        // This can be done independent of L1X level. All variables filled here are in groups
        // different from science_data.
        const vector<size_t> strt = {frm->l1x_iframe};
        const vector<size_t> cnt = {1};

        // Image attributes.
        netcdf_check(nc,var_image_time.putVar(strt,cnt,&frm->image_time));
        netcdf_check(nc,var_binning_table.putVar(strt,cnt,&binning_table_id));
        netcdf_check(nc,var_nr_coadditions.putVar(strt,cnt,&frm->nr_coadditions));
        netcdf_check(nc,var_exposure_time.putVar(strt,cnt,&frm->exposure_time));

    }

    return 0;

} // }}}

int L1X::write( // {{{
    Frame *frm, // Image to write L1X from.
    Settings_noise *set, // Noise settings (only for L1X with noise).
    CKD *ckd // Calibration key data.
)
{

    if (execute) {

        // Define start and count vectors, relevant for all the binned images.
        vector<size_t> strt = {frm->l1x_iframe,0};
        vector<size_t> cnt = {1,npix_binned};

        if (step == L1X_L1A) {
            // For the L1A step, we can no longer add noise, because
            // we are already in an integer format. The noise at L1A level
            // is done by hand in the program.cpp.

            // Write image.
            netcdf_check(nc,var_img.putVar(strt,cnt,frm->image_ints.data()));
        }

        // Noise estimates.
        vector<double> noise_img;
        vector<double> noise_spectra;
        vector<double> write_img;
        vector<double> write_spectra;
        // Noise on images.
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
            // It is not necessary to calculate the noise for raw and dark
            // if no noise realization is included. Currently, I am too
            // lazy to make the if-clause even more complicated.
            noise_img.resize(npix);
            for (size_t ipix=0 ; ipix<npix ; ipix++) {
                noise_img[ipix] = 1.0;
            }
            // Add noise realization.
            write_img = frm->image; // Copy the vector.
        }
        // Noise on spectra.
        if (
            step == L1X_FOV ||
            step == L1X_RAD
        ) {
            const size_t nel = ckd->dim_fov*DIM_POL*ckd->dim_detector_spec;
            noise_spectra.resize(nel);
            for (size_t iel=0 ; iel<nel ; iel++) {
                noise_spectra[iel] = 1.0;
            }
            // Add noise realization.
            write_spectra = frm->intens;
        }
        // Write images, relevant for all these steps.
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

            // Image and mask.
            if (
                step == L1X_RAW ||
                step == L1X_DARK ||
                step == L1X_NOISE ||
                step == L1X_NONLIN ||
                step == L1X_PRNU ||
                step == L1X_REBIN
            ) {

                // The image.
                netcdf_check(nc,var_img.putVar(strt,cnt,write_img.data()));
                // Mask.
                vector<uint8_t> writemask(npix);
                for (size_t ipix=0 ; ipix<npix ; ipix++) writemask[ipix] = ckd->mask[ipix]?1:0;
                netcdf_check(nc,var_img_mask.putVar(strt,cnt,writemask.data()));

            }
            // Image noise, only relevant from noise step onwards. And
            // also relevant for straylight, because the
            // binned noise image remains in these steps.
            if (
                step == L1X_NOISE ||
                step == L1X_NONLIN ||
                step == L1X_PRNU ||
                step == L1X_UNBIN ||
                step == L1X_STRAY ||
                step == L1X_REBIN
            ) {
                netcdf_check(nc,var_img_noise.putVar(strt,cnt,noise_img.data()));
            }

            // For dark, also write the image with dark current not yet subtracted.
            if (step == L1X_DARK) {
                netcdf_check(nc,var_img_with_current.putVar(strt,cnt,frm->image_with_current.data())); // Screw the noise here.
            }
        }

        // Write unbinned image and mask.
        if (
            step == L1X_UNBIN ||
            step == L1X_STRAY
        ) {
            // Ironically, we have the binning table for the unbinned
            // detector dimensions. It contains CKD and 2D detector
            // dimensions are always unbinned.
            vector<size_t> strt = {frm->l1x_iframe,0,0};
            vector<size_t> cnt = {1,ckd->dim_detector_spat,ckd->dim_detector_spec};
            netcdf_check(nc,var_img.putVar(strt,cnt,write_img.data()));

            // The mask is created 1D, but written 2D. There is binning (yet).
            vector<uint8_t> writemask(npix);
            for (size_t ipix=0 ; ipix<npix ; ipix++) writemask[ipix] = ckd->mask[ipix]?1:0;
            netcdf_check(nc,var_img_mask.putVar(strt,cnt,writemask.data()));
        }

        // Spectra.
        if (
            step == L1X_FOV ||
            step == L1X_RAD
        ) {

            vector<uint32_t> dims_spec(ckd->dim_fov,ckd->dim_detector_spec);
            // The size of the spectra, one write instruction with different start and count vectors.
            netcdf_check(nc,var_spectra_sizes.putVar({frm->l1x_iframe,0},{1,ckd->dim_fov},dims_spec.data()));

            // Start and count vectors for the serious vectors.
            vector<size_t> strt = {frm->l1x_iframe,0,0,0};
            vector<size_t> cnt = {1,ckd->dim_fov,DIM_POL,ckd->dim_detector_spec};

            // Spectra, the signal.
            netcdf_check(nc,var_spectra.putVar(strt,cnt,write_spectra.data()));
            // Noise.
            netcdf_check(nc,var_spectra_noise.putVar(strt,cnt,noise_spectra.data()));
            // Mask. No issues about incomplete spectra, just the conversion from a boolean to an unsigned byte.
            vector<uint8_t> writemask(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec,0);
            netcdf_check(nc,var_spectra_mask.putVar(strt,cnt,writemask.data()));

            // Repeat joke for the 'target' spectra for RAD.
            if (step == L1X_RAD) {

                // Spectra, the signal.
                netcdf_check(nc,var_spectra_target.putVar(strt,cnt,write_spectra.data()));
                // Noise.
                netcdf_check(nc,var_spectra_target_noise.putVar(strt,cnt,noise_spectra.data()));
                // Mask. No issues about incomplete spectra, just the conversion from a boolean to an unsigned byte.

                // Recycle the writemask.
                netcdf_check(nc,var_spectra_target_mask.putVar(strt,cnt,writemask.data()));
            }
        }
    }

    return 0;

} // }}}

