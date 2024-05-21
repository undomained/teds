// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "matrix.h"
#include "fourier.h"
#include "lininv.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "settings_proc.h"
#include "ckd.h"
#include "l1a_file_metadata.h"
#include "binningtable.h"
#include "l1x.h"
#include "l1a.h"

#include "cubic_spline.h"

namespace tango {

// Constructor.
L1A::L1A( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
// Destructor.
L1A::~L1A() {}

// Reads the metadata.
// This is low-data-load information that either.
// - is required to detemine the batching.
// - may in the future be required for batching.
// - is derived from data with a different time axis.
// - is so similar to other metadata that it is logical that they are read in the same routine.
// At least the detector images themselves will not be included here.
// However, besides batching, no processing is done with L1A instances of which only the
// metadata is read.
int L1A::read_metadata( // {{{
)
{

    // Set shared pointer to the NetCDF object.
    // The member file_meta is going to die at the end of the init routine of the L1A manager.
    nc = file_meta->nc; // Here on the L1A structure, it is safe for as long as we need it.
    filename = file_meta->filename; // For NetCDF-unrelated error handling.
    file_id = file_meta->file_id; // File identifier.

    // Re-open the fie.
    nc->open(nc->filename,NcFile::read);

    // Raw NetCDF contents, encapsulated by the relevant frame indices.
    vector<double> time_raw(nframe_inp);
    vector<double> exposure_time_raw(nframe_inp);
    vector<int> digital_offset_raw(nframe_inp);
    vector<int> nr_coadditions_raw(nframe_inp);
    vector<uint8_t> binning_table_raw(nframe_inp);
    // Group image_attributes.
    NcGroup grp;
    NcGroup grpnd;
    const vector<size_t> strt = {frame_offset};
    const vector<size_t> cnt = {nframe_inp};
    netcdf_check(nc,grp = nc->ncid->getGroup("image_attributes"));
    netcdf_check(nc,grp.getVar("image_time").getVar(strt,cnt,time_raw.data()));
    netcdf_check(nc,grp.getVar("exposure_time").getVar(strt,cnt,exposure_time_raw.data()));
    netcdf_check(nc,grp.getVar("nr_coadditions").getVar(strt,cnt,nr_coadditions_raw.data()));
    netcdf_check(nc,grp.getVar("binning_table").getVar(strt,cnt,binning_table_raw.data()));

    // Close the file again.
    nc->close();

    // The binning table should be consistent within an image.
    binning_table_id = binning_table_raw[0];
    for (size_t iframe=1 ; iframe<nframe_inp ; iframe++) {
        check_error(binning_table_raw[iframe] != binning_table_id,"Error: Inconsistent binning table within frames that belong to the same image.");
    }

    // Co-add the frames to one image.
    // Initialize aggregates to zero.
    exposure_time = 0.0;
    nr_coadditions = 0;
    for (size_t iframe=0 ; iframe<nframe_inp ; iframe++) {
        // Co-add the metadata.
        exposure_time += nr_coadditions_raw[iframe]*exposure_time_raw[iframe];
        nr_coadditions += nr_coadditions_raw[iframe];
    }
    // Exposure time should be averaged.
    exposure_time /= nr_coadditions;

    // We will always calculate the time that belongs to the entire image.

    // The time we will use is the average of the raw time for the first and last frame.
    // For flight data, where navigation data is relevant, the first and last frames are the same.
    // We will only use this time for navigation data, but we will save the time even if no
    // navigation data is requested.
    time = 0.5*(time_raw[0] + time_raw[nframe_inp-1]);

    act_angle = file_meta->act_angle;
    alt_angle = file_meta->alt_angle;

    // Copy L1X iformation.
    il1x_start = file_meta->il1x_start;
    l1x_dim_detector_spec = file_meta->l1x_dim_detector_spec;
    l1x_dim_detector_spat = file_meta->l1x_dim_detector_spat;
    l1x_dim_fov = file_meta->l1x_dim_fov;

    return 0;

} // }}}

// Reads the image.
// This is high-data-load.
// These are, thus, the detector image frames and eventually co-adds them.
int L1A::read_image( // {{{
    size_t a_ipix_start,
    size_t a_ipix_end
)
{

    writelog(log_verbose,"Reading image from '%s'.",filename.c_str());

    // Re-open the fie.
    nc->open(nc->filename,NcFile::read);

    // Copy chunk information.
    ipix_start = a_ipix_start;
    ipix_end = a_ipix_end;

    check_error(ipix_end < ipix_start,"Error: Negative image fraction selected, from %zu to %zu.",ipix_start,ipix_end);
    check_error(ipix_end > bin->npix,"Error: Image too small for selected fraction. Image size is %zu, selected fraction goes to %zu.",bin->npix,ipix_end);

    // Size of the fraction.
    size_t npix_frac = ipix_end-ipix_start;

    // Group: science data.
    NcGroup grp;
    netcdf_check(nc,grp = nc->ncid->getGroup("science_data"));

    // Allocate and read (binned) image and pixel mask.
    if (
        il1x_start == L1X_L1A ||
        il1x_start == L1X_RAW ||
        il1x_start == L1X_DARK ||
        il1x_start == L1X_NOISE ||
        il1x_start == L1X_NONLIN ||
        il1x_start == L1X_PRNU ||
        il1x_start == L1X_REBIN
    ) {

        // Ultimate image.
        image_read.resize(npix_frac);
        // Set pointer so that normal pixel index can be used.
        image = image_read.data() - ipix_start;

        // Custom pixel mask.
        pixelmask.resize(bin->npix,true);
        // Always fully defined, because bool* and vector<bool> do not have the desired relationship.
        // Pixels are dead unless they are proven to be alive.

        // Read the image. This is difficult for L1A, but easy for the
        // other ones. At least there is a detector images variable.

        if (il1x_start == L1X_L1A) {

            // Start and count for a binned image.
            const vector<size_t> strt = {frame_offset,ipix_start};
            const vector<size_t> cnt = {nframe_inp,npix_frac};

            // Read normal L1A if there is no L1X maturity.
            vector<int> raw(nframe_inp*npix_frac);
            netcdf_check(nc,grp.getVar("detector_images").getVar(strt,cnt,raw.data()));

            // Ultimate image.
            image_read.resize(npix_frac);
            // Set pointer so that normal pixel index can be used.
            image = image_read.data() - ipix_start;

            // Raw frames. Each pointer can be used with normal pixel index.
            vector<int *> rawframes(nframe_inp);
            for (size_t iframe=0 ; iframe<nframe_inp ; iframe++) {
                rawframes[iframe] = &raw[iframe*npix_frac] - ipix_start;
            }

            // Pointer to the image, which is frame zero.
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                image[ipix] = 0.0;
            }
            for (size_t iframe=0 ; iframe<nframe_inp ; iframe++) {
                for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                    image[ipix] += rawframes[iframe][ipix];
                }
            }
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                size_t binsize = bin->getBinsize(ipix);
                if (binsize == 0) image[ipix] = NC_FILL_DOUBLE;
                else image[ipix] /= nr_coadditions*binsize;
            }

            // Pixels inside the fraction are alive.
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                pixelmask[ipix] = false;
            }
        } else {

            // Start and count for a binned image.
            const vector<size_t> strt = {frame_offset,ipix_start};
            const vector<size_t> cnt = {nframe_inp,npix_frac}; // The nframe_inp will always be one for any L1X that is not L1A.
            // L1X maturity, at least there is an image in floating points.
            netcdf_check(nc,grp.getVar("detector_images").getVar(strt,cnt,image_read.data()));
            // Pixel mask is now also a variable.
            vector<uint8_t> maskread(npix_frac);
            netcdf_check(nc,grp.getVar("detector_images_mask").getVar(strt,cnt,maskread.data()));
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) pixelmask[ipix] = (maskread[ipix-ipix_start] != 0);
        }
    }
    // Read binned noise image.
    if (
        il1x_start == L1X_NOISE ||
        il1x_start == L1X_NONLIN ||
        il1x_start == L1X_PRNU ||
        il1x_start == L1X_UNBIN ||
        il1x_start == L1X_STRAY ||
        il1x_start == L1X_REBIN
    ) {
        // Noise. Similar to the image.
        // Start and count for a binned image.
        const vector<size_t> strt = {frame_offset,ipix_start};
        const vector<size_t> cnt = {nframe_inp,npix_frac}; // The nframe_inp will always be one for any L1X that is not L1A.
        noise_read.resize(npix_frac);
        noise = noise_read.data() - ipix_start;
        netcdf_check(nc,grp.getVar("detector_images_noise").getVar(strt,cnt,noise_read.data()));
    }
    // Read image with current if it exists. Only for L1X_DARK.
    if (il1x_start == L1X_DARK) {
        // Start and count for a binned image.
        const vector<size_t> strt = {frame_offset,ipix_start};
        const vector<size_t> cnt = {nframe_inp,npix_frac}; // The nframe_inp will always be one for any L1X that is not L1A.
        image_with_dark_current_read.resize(npix_frac);
        image_with_dark_current = image_with_dark_current_read.data() - ipix_start;
        netcdf_check(nc,grp.getVar("detector_images_with_current").getVar(strt,cnt,image_with_dark_current_read.data()));
    }

    // Unbinned image and mask, also stored in the same variables.
    // From here on, it is impossible to use a fractional L1A file.
    if (
        il1x_start == L1X_UNBIN ||
        il1x_start == L1X_STRAY
    ) {

        size_t npix_unbinned = l1x_dim_detector_spat * l1x_dim_detector_spec;
        const vector<size_t> strt = {frame_offset,0,0};
        const vector<size_t> cnt = {nframe_inp,l1x_dim_detector_spat,l1x_dim_detector_spec}; // The nframe_inp will always be one for any L1X that is not L1A.
        image_read.resize(npix_unbinned);
        image = image_read.data(); // No fractionation.
        netcdf_check(nc,grp.getVar("detector_images_2d").getVar(strt,cnt,image));
        pixelmask.resize(npix_unbinned);
        vector<uint8_t> maskread(npix_unbinned);
        netcdf_check(nc,grp.getVar("detector_images_2d_mask").getVar(strt,cnt,maskread.data()));
        for (size_t ipix=0 ; ipix<npix_unbinned ; ipix++) pixelmask[ipix] = (maskread[ipix] != 0);

    }

    // Spectra.
    if (
        il1x_start == L1X_FOV ||
        il1x_start == L1X_RAD
    ) {
        // Allocate spectra stash.
        l1x_stashed_spectra.resize(l1x_dim_fov);
        // Get spectrum size.
        NcVar var_spectra_sizes;
        netcdf_check(nc,var_spectra_sizes = grp.getVar("spectra_sizes"));
        NcVar var_spectra;
        netcdf_check(nc,var_spectra = grp.getVar(il1x_start == L1X_RAD?"spectra_target":"spectra"));
        NcVar var_spectra_noise;
        netcdf_check(nc,var_spectra_noise = grp.getVar(il1x_start == L1X_RAD?"spectra_target_noise":"spectra_noise"));
        NcVar var_spectra_mask;
        netcdf_check(nc,var_spectra_mask = grp.getVar(il1x_start == L1X_RAD?"spectra_target_mask":"spectra_mask"));
        vector<size_t> strt = {frame_offset,NC_FILL_UINT64,0,0};
        vector<size_t> cnt = {nframe_inp,1,DIM_POL,NC_FILL_UINT64};
        for (size_t ifov=0 ; ifov<l1x_dim_fov ; ifov++) {
            // Get the spectrum size (with inline start and count vectors).
            uint32_t dim_spec;
            netcdf_check(nc,var_spectra_sizes.getVar({frame_offset,ifov},{nframe_inp,1},&dim_spec)); // The nframe_inp is also a one.
            strt[1] = ifov;
            cnt[3] = dim_spec;
            Spectra &specs = l1x_stashed_spectra[ifov];
            specs.dim = dim_spec;
            specs.signal.resize(DIM_POL*dim_spec);
            specs.noise.resize(DIM_POL*dim_spec);
            specs.mask.resize(DIM_POL*dim_spec);
            vector<uint8_t> spectra_maskread(DIM_POL*dim_spec);
            netcdf_check(nc,var_spectra.getVar(strt,cnt,specs.signal.data()));
            netcdf_check(nc,var_spectra_noise.getVar(strt,cnt,specs.noise.data()));
            netcdf_check(nc,var_spectra_mask.getVar(strt,cnt,spectra_maskread.data()));
            for (size_t iel=0 ; iel<DIM_POL*dim_spec ; iel++) {
                specs.mask[iel] = (spectra_maskread[iel] != 0);
            }
        }
    }

    // Close the file again.
    nc->close();

    return 0;

} // }}}

int L1A::remove_image( // {{{
)
{

    writelog(log_verbose,"Removing image from '%s'.",filename.c_str());
    // This routines removes everything that is created in the read routine.
    // Even the insignificant pointer arrays, just for symmetry. The image
    // and its noise (and the pixel mask) are the most important ones.
    ipix_start = 0;
    ipix_end = 0;
    image_read = vector<double>();
    image = NULL;
    noise_read = vector<double>();
    noise = NULL;
    pixelmask = vector<bool>();
    // Image with dark current should never survive long. If it is read as
    // L1X, the calibration of the detector should remove the array
    // automatically.
    // Spectra stash for L1X input. These are not cleaned automatically.
    l1x_stashed_spectra = vector<Spectra>();

    return 0;

} // }}}

// Smooth over bad values. This is necessary for stray light which
// uses all the image pixels.
static auto fillHoles(const std::vector<bool>& pixelmask,
                      const int npix,
                      double* image) -> void
{
    for (int i {}; i < npix; ++i) {
        if (pixelmask[i] || std::abs(image[i]) > 1e8) {
            // Unless we are at either end of the image array, take
            // the average of the neighboring values.
            if (i == 0) {
                image[0] = image[1];
            } else if (i == npix - 1) {
                image[npix - 1] = image[npix - 2];
            } else {
                image[i] = (image[i - 1] + image[i + 1]) / 2;
            }
            // If there are two or more bad pixels side by side we
            // could still end up with a bad value after
            // averaging. Then we just set it to zero. To test it we
            // use a large value that is close to the NetCDF float
            // fill value.
            constexpr double bad_value { 1e30 };
            if (image[i] > bad_value) {
                image[i] = 0.0;
            }
        }
    }
}

// Calibrates a L1A file up to a certain level.
// For this level the CKD must have that level as well.
int L1A::calibrate_detector( // {{{
    CKD *ckd_gen, // Regular calibration key data (for straylight). Electronic detector steps are done by bin->binned_ckd, which may be the same CKD.
    Calibration_options &opt // Detector calibration options.
)
{
    auto time_ini = chrono::high_resolution_clock::now(); // Speed check

    // Start using the CKD from the binning table.
    CKD *ckd = bin->binned_ckd;

    // Often-used variable.
    size_t npix_frac = ipix_end-ipix_start;

    // For L1X images, assert the correct dimensions.
    if (
        il1x_start == L1X_UNBIN ||
        il1x_start == L1X_STRAY
    ) {
        // Unbinned detector shape.
        check_error(l1x_dim_detector_spat != ckd->dim_detector_spat,"Error: Spatial detector dimension unequal between CKD (%zu) and L1X input (%zu)",ckd->dim_detector_spat,l1x_dim_detector_spat);
        check_error(l1x_dim_detector_spec != ckd->dim_detector_spec,"Error: Spectral detector dimension unequal between CKD (%zu) and L1X input (%zu)",ckd->dim_detector_spec,l1x_dim_detector_spec);
    }

    // CKD pixel mask. It is created during the dim step.
    if (il1x_start < L1X_RAW && ckd->lev > LEVEL_DIMCAL) {
        // What happens here is that dead pixels infect the entire bin.
        // And during unbinning, they will all be ill. The infection is
        // resolved during creation of the binned CKD.
        for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
            pixelmask[ipix] = ckd->mask[ipix];
        }
        // Write L1X output.
        handle(l1x[L1X_RAW]->write(this));
        // Detector dimensions are checked when reading the binning table.
    }

    // Dark calibartion.
    if (il1x_start < L1X_DARK && ckd->lev > LEVEL_DARKCAL) {
        // Save mage with the dark current, so that the noise can be calculated with that one.
        image_with_dark_current_read.resize(npix_frac); // This will be subtracted with a delay. Keep it zero if dark is skipped.
        image_with_dark_current = image_with_dark_current_read.data() - ipix_start; // Note that this pointer is a member, but it becomes dangling after calibrate_detector. Besides this routine, the pointer is only used in the L1X.
        if (!ckd->dark_skip && opt.dark_apply) {

            // Evaluate the fit. The B-splines are fixed to a polynomial where the first term
            // is the value at nominal temperature.
            vector<double> knots = {ckd->dark_nominal_temperature,ckd->dark_nominal_temperature+1.0};
            Bspline b(this,ckd->dim_dark_order,2,knots.data());
            vector<double> terms(ckd->dim_dark_order);

            double temp_detector { 300.0 };
            handle(b.jaccalc(1,&temp_detector,terms.data()));
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                if (!pixelmask[ipix]) {
                    for (size_t iorder=0 ; iorder<ckd->dim_dark_order ; iorder++) {
                        // Subtract dark offset.
                        image[ipix] -= ckd->dark_offset[ipix*ckd->dim_dark_order+iorder]*terms[iorder];
                    }
                    // Make copy with dark current.
                    image_with_dark_current[ipix] = image[ipix];
                    // Update real image if dark current must be subtracted.
                    // This is false for the noisecal calibration generation,
                    // so they can use the regular image. They should, because
                    // the pointer image_with_dark_current is not available
                    // in noisecal.
                    // For application of the noise, the image_with_dark_current
                    // is available, and also for the L1X.
                    for (size_t iorder=0 ; iorder<ckd->dim_dark_order ; iorder++) {
                        image[ipix] -= ckd->dark_current[ipix*ckd->dim_dark_order+iorder]*terms[iorder]*exposure_time;
                    }
                }
            }
        } else {
            image_with_dark_current = image;
            std::cout << "Skipping dark offset and current calibration"
                      << std::endl;
        }
        handle(l1x[L1X_DARK]->write(this));
    }

    // Noise calibration.
    if (il1x_start < L1X_NOISE && ckd->lev > LEVEL_NOISECAL) {
        noise_read.resize(npix_frac);
        noise = noise_read.data() - ipix_start;
        if (!ckd->noise_skip) {
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                if (pixelmask[ipix]) noise[ipix] = NC_FILL_DOUBLE;
                else {
                    double variance = ckd->noise_g[ipix]*image_with_dark_current[ipix] + pow(ckd->noise_n[ipix],2.0);
                    if (variance <= 0.0) {
                        image[ipix] = NC_FILL_DOUBLE;
                        noise[ipix] = NC_FILL_DOUBLE;
                        pixelmask[ipix] = true;
                    } else {
                        noise[ipix] = sqrt(variance) / sqrt(nr_coadditions*bin->getBinsize(ipix));
                    }
                }
            }
        } else {
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                if (pixelmask[ipix]) noise[ipix] = NC_FILL_DOUBLE;
                else noise[ipix] = 1.0;
            }
        }
        handle(l1x[L1X_NOISE]->write(this));
    }

    // The image without dark current correction is no longer relevant
    // fro here onwards. We throw away the pointer because it should
    // no longer be used.
    image_with_dark_current_read = vector<double>(); // Saves memory.
    image_with_dark_current = NULL;

    // Non-linearity correction.
    if (il1x_start < L1X_NONLIN && ckd->lev > LEVEL_NONLINCAL && opt.nonlin_apply) { // {{{
        if (!ckd->nonlin_skip) {
            // Choose the relevant exposure time.
            double mindiff = abs(ckd->nonlin_exptimes[0]-exposure_time);
            size_t iexptime_best = 0;
            for (size_t iexptime=1 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
                if (abs(ckd->nonlin_exptimes[iexptime]-exposure_time) < mindiff) {
                    mindiff = abs(ckd->nonlin_exptimes[iexptime]-exposure_time);
                    iexptime_best = iexptime;
                }
            }
            // Use the B-spline to evaluate the fit on the desired input (the signal).
            Bspline b(this,ckd->nonlin_order,ckd->dim_nonlin_knot,ckd->nonlin_knots.data());
            vector<size_t> ideriv = {0,1};
            double *nonlin_fit_cur = &ckd->nonlin_fit[ipix_start*ckd->dim_nonlin_exptime*ckd->dim_nonlin_spline+iexptime_best*ckd->dim_nonlin_spline];
            double *nonlin_signal_scale_cur = &ckd->nonlin_signal_scale[iexptime_best];
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                if (!pixelmask[ipix]) {
                    double &orig = image[ipix];
                    double signal_cur = orig;
                    bool converged = false;
                    for (size_t iter=0 ; iter<opt.nonlin_niter ; iter++) {
                        vector<double> model(2);
                        vector<double> abscissa(2);
                        abscissa[0] = signal_cur / (*nonlin_signal_scale_cur);
                        abscissa[1] = abscissa[0];
                        handle(b.jacapply(2,abscissa.data(),ideriv.data(),nonlin_fit_cur,model.data()));
                        double gain = *nonlin_signal_scale_cur / model[1];
                        double addition = (orig - model[0]) * gain;
                        signal_cur += addition;
                        if (abs(addition) < opt.nonlin_tol) {
                            orig = signal_cur;
                            noise[ipix] *= abs(gain);
                            converged = true;
                            break;
                        }
                    }
                    // Flag pixel if non-linearity inversion did not converge.
                    if (!converged) {
                        orig = NC_FILL_DOUBLE;
                        noise[ipix] = NC_FILL_DOUBLE;
                        pixelmask[ipix] = true;
                    }
                }
                nonlin_fit_cur += ckd->dim_nonlin_spline*ckd->dim_nonlin_exptime;
                nonlin_signal_scale_cur += ckd->dim_nonlin_exptime;
            }
        }
        handle(l1x[L1X_NONLIN]->write(this));
    } // }}}

    // Pixel-response non-uniformity.
    if (il1x_start < L1X_PRNU && ckd->lev > LEVEL_PRNUCAL) {
        if (!ckd->prnu_skip && opt.prnu_apply) {
            // Correct for pixel response non-uniformity.
            for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                if (!pixelmask[ipix]) {
                    image[ipix] /= ckd->prnu_prnu[ipix];
                    noise[ipix] /= ckd->prnu_prnu[ipix];
                }
            }
        } else {
            std::cout << "Skipping PRNU calibration" << std::endl;
        }
        handle(l1x[L1X_PRNU]->write(this));
    }
    // Straylight step. This consists of three substeps. Unbinning,
    // straylight correction, rebinning. Several special cases can apply/
    // It is possible that straylight is skipped or not selected for
    // the specific case, because it is a qualitatieve measurement. In such
    // a case, the sequence can be skipped unless the L1X input is somewhere
    // in the middle, e.g. L1X_UNBIN. Then, the data should be properly
    // rebinned. If the CKD does not have straylight yet, unbinning is
    // useless and should not be considered. However, if you only want to
    // write L1X output, what to do? I think it is best to make straylight
    // CKD, possibly with the skip flag and then write L1X output. Actually,
    // unbinning and rebinning should always be done, except that it is a
    // waste if there is no straylight correction. I do not know how much
    // this matters, because operationally, there will be straylight
    // correction. And for testing, how much time does this cost. It saves
    // a lot of logistics just to always perform the unbinning/rebinning.
    // TODO: Reconsider this whenever something changes in straylight,
    // e.g. if it is no longer necessary, or when unbinning an rebinning
    // appears to take so much time that it is irritating during testing.

    // Conclusion for now. CKD must contain straylight (eventually skip).
    // Then, all steps are done, although the straylight part is skipped
    // if straylight CKD says skip or the calibration options say skip
    // straylight.
    if (ckd->lev > LEVEL_STRAYCAL) {
        // Un-binning.
        // We will not unbin the noise. That is bullshit. Noise unbinned
        // becomes correlated misery. And straylight does not do anything
        // with noise (luckily). Good point, we will demolish unbinned noise
        // from the L1X, and we have to write binned noise, also for these
        // L1X.
        if (il1x_start < L1X_UNBIN) { // {{{
            if (!bin->trivial) {

                size_t &npix_unbinned = ckd_gen->npix; // Reference to unbinned number of pixels (or bin->npix_unbinned, is the same).

                // Apply the unbinning. This writes full-detector images.
                // If you use fractions, parts of the field will not be filled and
                // the pixels are masked. We will write from the first up to the
                // last detector pixel in the L1X file, so the gaps will be written
                // as fillvalue and the pixels before and after will be left alone
                // and stay fillvalue.

                // Apply binning table to convert it to the original detector.
                // TODO: Do this better. The placeholder just copies the value to each
                // of the involved pixels.
                vector<double> image_unbinned(npix_unbinned,NC_FILL_DOUBLE);
                vector<bool> pixelmask_unbinned(npix_unbinned,true); // Pixels are dead until they are flagged as alive.
                // By the way, pixels flagged as bad in the CKD will get their data, but are flagged as dead. Non-included pixels are set to fillvalue (and flagged as dead).
                for (size_t ipix_unbinned=0 ; ipix_unbinned<npix_unbinned ; ipix_unbinned++) {
                    uint32_t &ipix = bin->pixelpointer[ipix_unbinned];
                    if (ipix >= bin->npix) continue; // Thrown-away pixel.
                    if (ipix < ipix_start) continue; // Not yet in fraction (not recommended to use fractions anyway).
                    if (ipix >= ipix_end) {
                        break; // Passed your fraction. This is not a thrown-away pixel, because that is caught before.
                    }
                    image_unbinned[ipix_unbinned] = image[ipix];
                    pixelmask_unbinned[ipix_unbinned] = pixelmask[ipix];
                }
                // Overwrite the image and the pixel mask with the
                // full-image one.
                // This throws away the old vectors and puts the new information
                // in there.
                image_read = image_unbinned;
                pixelmask = pixelmask_unbinned;
                // Set pointer accordingly. They were lost while overwriting.
                image = image_read.data();
                // Everything is full-frame, so there is no ipix_start anymore.
            } else if (ipix_start != 0 || ipix_end != ckd->npix) {
                // Add fillvalues for non-included pixels.
                vector<double> image_unbinned(ckd->npix,NC_FILL_DOUBLE);
                vector<bool> pixelmask_unbinned(ckd->npix,true); // Pixels are dead until they are flagged as alive.
                for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {
                    image_unbinned[ipix] = image[ipix];
                    pixelmask_unbinned[ipix] = pixelmask[ipix];
                }
                // Overwrite the image and the pixel mask with the
                // full-image one.
                // This throws away the old vectors and puts the new information
                // in there.
                image_read = image_unbinned;
                pixelmask = pixelmask_unbinned;
                // Set pointer accordingly. They were lost while overwriting.
                image = image_read.data();
                // Everything is full-frame, so there is no ipix_start anymore.
            } // Else, do nothing.
            // Write L1X output (after unbinning).
            handle(l1x[L1X_UNBIN]->write(this));
        } // }}}
        // Straylight correction itself.
        if (il1x_start < L1X_STRAY) { // {{{
            CKD *ckd = ckd_gen; // This temporarily hides the binned CKD,
            if (ckd->stray_skip || opt.stray_van_cittert_steps == 0) {
                std::cout << "Skipping stray light correction" << std::endl;
            } else {
                fillHoles(pixelmask, static_cast<int>(ckd->npix), image);
            }
            if (!ckd->stray_skip) {
                // "Ideal" image, i.e. the one without stray light
                std::vector<double> image_ideal(ckd->npix);
                for (int i {}; i < ckd->npix; ++i) {
                    image_ideal[i] = image[i];
                }
                // Image multiplied by kernel weights
                std::vector<double> image_weighted(ckd->npix);
                // Result of taking a convolution
                std::vector<double> conv_result(ckd->npix);
                // Van Cittert algorithm
                for (int i_vc {}; i_vc < opt.stray_van_cittert_steps; ++i_vc) {
                    std::ranges::fill(conv_result, 0.0);
                    for (int i_kernel {}; i_kernel < ckd->stray.n_kernels;
                         ++i_kernel) {
                        // Each iteration starts with image_ideal as the best
                        // current estimate.
                        for (int i {}; i < ckd->npix; ++i) {
                            image_weighted[i] =
                              image_ideal[i] * ckd->stray.weights[i_kernel][i];
                        }
                        // Number of rows and column in this subimage
                        const int image_n_rows {
                            ckd->stray.edges[i_kernel * box::n + box::t]
                            - ckd->stray.edges[i_kernel * box::n + box::b]
                        };
                        const int image_n_cols {
                            ckd->stray.edges[i_kernel * box::n + box::r]
                            - ckd->stray.edges[i_kernel * box::n + box::l]
                        };
                        std::vector<double> sub_image(image_n_rows
                                                      * image_n_cols);
                        for (int i {}; i < image_n_rows; ++i) {
                            for (int j {}; j < image_n_cols; ++j) {
                                sub_image[i * image_n_cols + j] =
                                  image_weighted
                                  [(i + ckd->stray.edges[i_kernel * box::n
                                                         + box::b])
                                   * ckd->stray.n_spectral
                                   + j + ckd->stray.edges[i_kernel * box::n
                                                          + box::l]];
                            }
                        }
                        // Result of taking a convolution using one of the
                        // subimages and kernels
                        std::vector<double> conv_result_sub(ckd->npix);
                        convolve_fft(image_n_rows,
                                     image_n_cols,
                                     sub_image,
                                     ckd->stray.kernel_rows[i_kernel],
                                     ckd->stray.kernel_cols[i_kernel],
                                     ckd->stray.kernel_fft_sizes[i_kernel],
                                     ckd->stray.kernels_fft[i_kernel],
                                     conv_result_sub);
                        // The full convolution is a sum over all convolutions
                        for (int i {}; i < image_n_rows; ++i) {
                            for (int j {}; j < image_n_cols; ++j) {
                                conv_result
                                  [(i + ckd->stray.edges[i_kernel * box::n
                                                         + box::b])
                                   * ckd->stray.n_spectral
                                   + j + ckd->stray.edges[i_kernel * box::n
                                                          + box::l]] +=
                                  conv_result_sub[i * image_n_cols + j];
                            }
                        }
                    }
                    for (int i {}; i < ckd->npix; ++i) {
                        image_ideal[i] = (image[i] - conv_result[i])
                                         / (1 - ckd->stray.eta[i]);
                    }
                }
                for (int i {}; i < ckd->npix; ++i) {
                    image[i] = image_ideal[i];
                }
            }
            handle(l1x[L1X_STRAY]->write(this));
        } // }}}
        // Re-binning.
        if (il1x_start < L1X_REBIN) { // {{{
            // It seems a pity that we rebin a pixel mask with which nothing
            // has happened. This can be optimized, but that needs us to
            // save an additional mask and possibly, the unbinned mask is
            // read from the L1X. In the first implementation, we avoid the
            // logistic minefield and just do the binning and rebinning.
            // Optimization is foreseen if this step appears to take a
            // significant amount of time.

            // Fractional L1A should not work anymore here.
            // So, only the unbinning has to be performed.
            if (!bin->trivial) {
                // The image should be averaged. The mask will be true if
                // ay sub-pixel is dead. If straylight correction does not
                // kill off pixels, this means that the entire bin is alive
                // or dead.
                vector<double> image_rebinned(ckd->npix,0.0);
                vector<bool> pixelmask_rebinned(ckd->npix,false);
                for (size_t ipix_unbinned=0 ; ipix_unbinned<ckd_gen->npix ; ipix_unbinned++) {
                    uint32_t &ipix = bin->pixelpointer[ipix_unbinned];
                    if (ipix <= ckd->npix) {
                        image_rebinned[ipix] += image[ipix_unbinned];
                        if (pixelmask[ipix_unbinned]) pixelmask_rebinned[ipix] = true;
                    }
                }
                for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
                    image_rebinned[ipix] /= bin->getBinsize(ipix);
                }

                // Overwrite the full-image image and pixel mask with
                // the rebinned ones.
                image_read = image_rebinned;
                pixelmask = pixelmask_rebinned;

                // Restore pointer.
                image = image_read.data();
            }
            handle(l1x[L1X_REBIN]->write(this));
        } // }}}
    }
    auto time_fin = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = time_fin-time_ini;
    writelog(log_verbose,"Time for detector calibration: %.5f s.",elapsed.count());
    return 0;

} // }}}

// Extracts S+ and S- spectrum from detector image.
int L1A::extract( // {{{
    size_t ifov, // Field-of-view index.
    const Calibration_options& opt,
    Spectra &specs // Output spectra.
)
{

    CKD *ckd = bin->binned_ckd;

    // For L1X input, check dimension consistency.
    if (
        il1x_start == L1X_FOV ||
        il1x_start == L1X_RAD
    ) {
        check_error(l1x_dim_fov != ckd->dim_fov,"Error: Number of spatial samples per image is unqual between CKD (%zu) and L1X input (%zu)",ckd->dim_fov,l1x_dim_fov);
        check_error(l1x_stashed_spectra[ifov].dim != ckd->fov_dims_spec[ifov],"Error: Number of spectral samples of extracted spectrum of FOV number %zu is unequal between CKD (%zu) and L1X input (%zu)",ifov,ckd->fov_dims_spec[ifov],l1x_stashed_spectra[ifov].dim);
    }

    // FOV-calibration.
    if (ckd->lev > LEVEL_FOVCAL) {

        // L1X maturity lower than FOV, this step is new and should be executed
        // like normal.
        if (il1x_start < L1X_FOV) {

            // Determine spectra size.
            specs.dim = ckd->fov_dims_spec[ifov];

            // Allocate the spectra.
            specs.signal.resize(DIM_POL*specs.dim);
            specs.noise.resize(DIM_POL*specs.dim);
            specs.mask.resize(DIM_POL*specs.dim);

            // Offset pointers over DIM_POL*fov_dims_spec elements.
            size_t *fov_ipix1_cur = &ckd->fov_ipix1[ckd->fov_iel_start[ifov]];
            size_t *fov_ipix2_cur = &ckd->fov_ipix2[ckd->fov_iel_start[ifov]];
            double *fov_weight1_cur = &ckd->fov_weight1[ckd->fov_iel_start[ifov]];

            for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                for (size_t ispec=0 ; ispec<specs.dim ; ispec++) {
                    size_t &ipix_left = fov_ipix1_cur[ipol*specs.dim+ispec];
                    size_t &ipix_right = fov_ipix2_cur[ipol*specs.dim+ispec];
                    double &weightleft = fov_weight1_cur[ipol*specs.dim+ispec];
                    if (ipix_left >= ckd->npix && ipix_right >= ckd->npix) {
                        specs.signal[ipol*specs.dim+ispec] = NC_FILL_DOUBLE;
                        specs.noise[ipol*specs.dim+ispec] = NC_FILL_DOUBLE;
                        specs.mask[ipol*specs.dim+ispec] = true;
                    } else {
                        if (pixelmask[ipix_left]) {
                            if (pixelmask[ipix_right]) {
                                specs.signal[ipol*specs.dim+ispec] = NC_FILL_DOUBLE;
                                specs.noise[ipol*specs.dim+ispec] = NC_FILL_DOUBLE;
                                specs.mask[ipol*specs.dim+ispec] = true;
                            } else {
                                specs.signal[ipol*specs.dim+ispec] = image[ipix_right];
                                specs.noise[ipol*specs.dim+ispec] = noise[ipix_right];
                            }
                        } else {
                            if (pixelmask[ipix_right]) {
                                specs.signal[ipol*specs.dim+ispec] = image[ipix_left];
                                specs.noise[ipol*specs.dim+ispec] = noise[ipix_left];
                            } else {
                                specs.signal[ipol*specs.dim+ispec] = weightleft * image[ipix_left] + (1.0-weightleft) * image[ipix_right];
                                specs.noise[ipol*specs.dim+ispec] = weightleft * noise[ipix_left] + (1.0-weightleft) * noise[ipix_right];
                            }
                        }
                    }
                }
            }

            // Aim pointer to spectra for L1X.
            l1x_specs = &specs;
            handle(l1x[L1X_FOV]->write(this,ifov));
        }
        // L1X maturity is at FOV. Look up the spectra from the L1X stash.
        if (il1x_start == L1X_FOV) {
            specs = l1x_stashed_spectra[ifov]; // Copies all the contents.
        }
    }

    if (ckd->lev > LEVEL_RADCAL) {

        // If L1X maturity is not yet radiometrically calibrated, do what is normal.
        if (il1x_start < L1X_RAD) {

            if (!ckd->rad_skip && opt.rad_apply) {
                // We multiply by the calibration factor and divide by the exposure time, to make it radiance.
                double *rad_cur = &ckd->rad_spectra[ckd->fov_iel_start[ifov]];
                for (size_t iel=0 ; iel<DIM_POL*specs.dim ; iel++) {
                    if (specs.mask[iel]) continue;
                    if (rad_cur[iel] == NC_FILL_DOUBLE) {
                        specs.signal[iel] = NC_FILL_DOUBLE;
                        specs.noise[iel] = NC_FILL_DOUBLE;
                        specs.mask[iel] = true;
                    } else {
                        specs.signal[iel] *= rad_cur[iel] / exposure_time;
                        specs.noise[iel] *= rad_cur[iel] / exposure_time;
                    }
                }
            } else if (ifov == 0) {
                std::cout << "Skipping radiometric calibration" << std::endl;
            }

            // Now, interpolate on the common wavelength grid.
            double *wave_cur = &ckd->wave_spectra[ckd->fov_iel_start[ifov]];
            double *wave_target_cur = &ckd->wave_target[ckd->fov_iel_start[ifov]/DIM_POL];
            for (size_t ispec=0 ; ispec<specs.dim ; ispec++) {
                wave_target_cur[ispec] = wave_cur[ispec];
            }
            // Copy original spectra. This means that specs_cpy remains
            // how the spectra where before interpolation. And specs ends
            // up as the interpolated spectra.
            Spectra specs_cpy = specs;
            // Aim pointer to spectra for L1X.
            l1x_specs = &specs_cpy;
            l1x_specs_target = &specs;
            handle(l1x[L1X_RAD]->write(this,ifov));

        }
        // If L1X maturity is the radiometrically calibrated spectra, read the stashed spectra.
        // The read routine already carefully read the 'target' spectra.
        if (il1x_start == L1X_RAD) {
            specs = l1x_stashed_spectra[ifov]; // Copy all the contents.
        }

    }

    // Emphasize that the L1X pointers should no longer be used.
    l1x_specs = NULL;
    l1x_specs_target = NULL;

    return 0;

} // }}}

} // namespace tango
