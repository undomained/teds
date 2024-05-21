#include "header.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"

namespace tango {

// Constructor only passing argument to the parent.
CKD::CKD( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
CKD::~CKD() {}

// Read the CKD file and fills the structure as far as the CKD is there.
int CKD::read( // {{{
    Settings_main &set,
    level_t lev_target,
    bool write
)
{

    // Start with the first level. The CKD is such that you should start
    // running dimcal.
    lev = LEVEL_DIMCAL;

    // All CKD of earlier steps are read and meanwhile there are written in the
    // output CKD. There is no output CKD if the earliest step is L1B or L1C.
    // There is no input CKD if the earliest step is DIM.

    // The input NetCDF will be opened and closed in this module. The output
    // NetCDF will be opened and remained open and the objects are stored on the
    // CKD structure.
    NetCDF_object nc_in(this); // Will be opened if level > LEVEL_DIMCAL.

    // While reading any step, the contents are automatically written to the
    // output NetCDF. Also, unused (diagnostic) CKD is read and written, but not
    // stored on the CKD structure.

    // Create output file if relevant.
    write_ckd = write; // Copy argument to member variable.
    if (write_ckd) {
        nc_ckd = make_unique<NetCDF_object>(this);
        handle(nc_ckd->open(set.ckd_file_out,NcFile::replace));
    }

    // Note that the level is the first level of CKD that will be added. So for
    // CKD to be already in the database, the first action level must be higher than
    // (not equal to) the relevant level.

    // Dimension: Get the dimensions, so that the CKD structure knows what it
    // should know (everything).
    if (lev_target > LEVEL_DIMCAL) { // {{{

        // Open input NetCDF file.
        handle(nc_in.open(set.ckd_file_in,NcFile::read));

        // Read dimension lengths of trivial dimensions and check they match the fixed ones.
        size_t dimcheck;
        netcdf_check(&nc_in,dimcheck = nc_in.ncid->getDim("vector_elements").getSize());
        check_error(dimcheck != DIM_VEC,"Error: Vector dimnension in CKD unequal to %d.",DIM_VEC);

        // Read non-trivial dimensions that were defined during dimensions calibration.
        netcdf_check(&nc_in,dim_detector_spat = nc_in.ncid->getDim("spatial_detector_pixels").getSize());
        netcdf_check(&nc_in,dim_detector_spec = nc_in.ncid->getDim("spectral_detector_pixels").getSize());
        netcdf_check(&nc_in,dim_vp = nc_in.ncid->getDim("number_of_views").getSize());
        // Number of pixels in general CKD is the product of length and width.
        // This is unbinned CKD.
        npix = dim_detector_spec*dim_detector_spat;

        // There is always a pixel mask, possibly not yet in its final state.
        // That is because the pixel mask is updated by all the detector calibration
        // steps.

        // Read the pixel mask.
        vector<uint8_t> maskread(npix); // NetCDF has no bools, we want bools.
        netcdf_check(&nc_in,nc_in.ncid->getVar("mask").getVar(maskread.data()));
        // Convert NetCDF byte to bool.
        mask.resize(npix);
        for (size_t ipix=0 ; ipix<npix ; ipix++) mask[ipix] = maskread[ipix] == 1;
        // Count number of living pixels (for convenience).
        nliving = 0;
        for (size_t ipix=0 ; ipix<npix ; ipix++) if (!mask[ipix]) nliving++;

        // Read the viewport mask. A viewport mask is there to turn off
        // viewports so that test calculations can be done more quickly.
        vector<uint8_t> vp_maskread(dim_vp); // NetCDF has no bools, we want bools.
        netcdf_check(&nc_in,nc_in.ncid->getVar("vp_mask").getVar(vp_maskread.data()));
        vp_mask.resize(dim_vp);
        for (size_t ivp=0 ; ivp<dim_vp ; ivp++) vp_mask[ivp] = vp_maskread[ivp] == 1;

    } // }}}
    if (lev_target > LEVEL_DARKCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("DARK"));

        // Define diagnostic variable.
        vector<double> dark_chi2;

        uint8_t intapply = set.dark_apply;
        dark_apply = intapply == 1;

        if (dark_apply) {

            // Read the dark correction.
            netcdf_check(&nc_in,dim_dark_order = grp.getDim("dark_number_of_coefficients").getSize());
            dark_offset.resize(npix*dim_dark_order);
            dark_current.resize(npix*dim_dark_order);
            netcdf_check(&nc_in,grp.getVar("dark_offset").getVar(dark_offset.data()));
            netcdf_check(&nc_in,grp.getVar("dark_current").getVar(dark_current.data()));
            netcdf_check(&nc_in,grp.getVar("dark_nominal_temperature").getVar(&dark_nominal_temperature));

            // Read useless stuff to copy it to the output.
            if (write_ckd) {
                // Shape temporary array.
                dark_chi2.resize(npix);
                // Read useless data.
                netcdf_check(&nc_in,grp.getVar("dark_chi2").getVar(dark_chi2.data()));
                // Set member pointer.
                diag_dark_chi2 = dark_chi2.data();
            }
        }
    } // }}}
    if (lev_target > LEVEL_NOISECAL) { // {{{

        // Noise CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("NOISE"));

        uint8_t intapply = set.noise_apply;
        noise_apply = intapply == 1;

        if (noise_apply) {

            // Shape the noise CKD.
            noise_g.resize(npix);
            noise_n.resize(npix);

            // Read the noise CKD.
            netcdf_check(&nc_in,grp.getVar("noise_g").getVar(noise_g.data()));
            netcdf_check(&nc_in,grp.getVar("noise_n").getVar(noise_n.data()));
        }
    } // }}}
    if (lev_target > LEVEL_NONLINCAL) { // {{{

        // Non-linearity CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("NON_LINEARITY"));

        // Define diagnostic variables.
        vector<double> nonlin_lin_slope;
        vector<double> nonlin_lin_chi2;
        vector<double> nonlin_chi2;

        uint8_t intapply = set.nonlin_apply;
        nonlin_apply = intapply == 1;

        if (nonlin_apply) {

            netcdf_check(&nc_in,dim_nonlin_exptime = grp.getDim("nonlin_number_of_exposure_times").getSize());
            netcdf_check(&nc_in,dim_nonlin_spline = grp.getDim("nonlin_number_of_coefficients").getSize());
            netcdf_check(&nc_in,dim_nonlin_knot = grp.getDim("nonlin_number_of_knots").getSize());

            // Shape CKD-owned arrays.
            nonlin_knots.resize(dim_nonlin_knot);
            nonlin_fit.resize(npix*dim_nonlin_exptime*dim_nonlin_spline);
            nonlin_exptimes.resize(dim_nonlin_exptime); // Exposure times for which non-linearity is calculated.
            nonlin_signal_scale.resize(npix*dim_nonlin_exptime); // Scaling signal for the abscissa of the non-linearity fit.

            int intorder; // For integer type conversion.
            netcdf_check(&nc_in,grp.getVar("nonlin_order").getVar(&intorder));
            nonlin_order = (size_t) intorder;
            netcdf_check(&nc_in,grp.getVar("nonlin_knots").getVar(nonlin_knots.data()));
            netcdf_check(&nc_in,grp.getVar("nonlin_exptimes").getVar(nonlin_exptimes.data()));
            netcdf_check(&nc_in,grp.getVar("nonlin_signal_scale").getVar(nonlin_signal_scale.data()));
            netcdf_check(&nc_in,grp.getVar("nonlin_fit").getVar(nonlin_fit.data()));

            // Read all useless stuff if it has to be copied to the output CKD file.
            if (write_ckd) {
                // Shape temporary arrays.
                nonlin_lin_slope.resize(npix);
                nonlin_lin_chi2.resize(npix);
                nonlin_chi2.resize(dim_nonlin_exptime*npix);
                // Read useless data.
                netcdf_check(&nc_in,grp.getVar("nonlin_lin_slope").getVar(nonlin_lin_slope.data()));
                netcdf_check(&nc_in,grp.getVar("nonlin_lin_chi2").getVar(nonlin_lin_chi2.data()));
                netcdf_check(&nc_in,grp.getVar("nonlin_chi2").getVar(nonlin_chi2.data()));
                // Set member pointers to temporary diagnostic data.
                diag_nonlin_lin_slope = nonlin_lin_slope.data();
                diag_nonlin_lin_chi2 = nonlin_lin_chi2.data();
                diag_nonlin_chi2 = nonlin_chi2.data();
            }
        }
    } // }}}
    if (lev_target > LEVEL_PRNUCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("PRNU"));

        uint8_t intapply = set.prnu_apply;
        prnu_apply = intapply == 1;

        if (prnu_apply) {
            // Pixel response non-uniformity.
            prnu_prnu.resize(npix);
            netcdf_check(&nc_in,grp.getVar("prnu_prnu").getVar(prnu_prnu.data()));
        }
    } // }}}
    if (lev_target > LEVEL_STRAYCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("STRAYLIGHT"));

        uint8_t intapply = set.stray_apply;
        stray_apply = intapply == 1;

        if (stray_apply) {
            stray.n_kernels = static_cast<int>(grp.getDim("kernels").getSize());
            stray.n_spatial = static_cast<int>(grp.getDim("spatial").getSize());
            stray.n_spectral =
              static_cast<int>(grp.getDim("spectral").getSize());
            stray.kernel_rows.resize(stray.n_kernels);
            grp.getVar("kernel_rows").getVar(stray.kernel_rows.data());
            stray.kernel_cols.resize(stray.n_kernels);
            grp.getVar("kernel_cols").getVar(stray.kernel_cols.data());
            stray.kernel_fft_sizes.resize(stray.n_kernels);
            grp.getVar("kernel_fft_sizes")
              .getVar(stray.kernel_fft_sizes.data());
            stray.kernels_fft.resize(stray.n_kernels);
            for (int i {}; i < stray.n_kernels; ++i) {
                stray.kernels_fft[i].resize(stray.kernel_fft_sizes[i]);
            }
            std::vector<size_t> start { 0 };
            for (int i_kernel {}; i_kernel < stray.n_kernels; ++i_kernel) {
                // Need to multiply by 2 because of complex numbers
                const std::vector<size_t> count { 2 * static_cast<size_t>(
                      stray.kernel_fft_sizes.at(i_kernel))
                };
                grp.getVar("kernels_fft")
                  .getVar(start,
                          count,
                          stray.kernels_fft.at(i_kernel).data());
                start.front() += count.front();
            }
            stray.eta.resize(stray.n_spatial * stray.n_spectral);
            grp.getVar("eta").getVar(stray.eta.data());
            stray.weights.resize(
              stray.n_kernels,
              std::vector<double>(stray.n_spatial * stray.n_spectral));
            std::vector<double> buf(stray.n_kernels * stray.n_spatial
                                    * stray.n_spectral);
            grp.getVar("weights").getVar(buf.data());
            for (int i {}; i < stray.n_kernels; ++i) {
                for (int j {}; j < stray.n_spatial * stray.n_spectral;
                     ++j) {
                    stray.weights[i][j] =
                      buf[i * stray.n_spatial * stray.n_spectral + j];
                }
            }
            stray.edges.resize(4 * stray.n_kernels);
            grp.getVar("edges").getVar(stray.edges.data());
        }
    } // }}}
    if (lev_target > LEVEL_FOVCAL) { // {{{

        // Field of view calibration.

        // This dimension is in the root, because later CKD steps
        // also use this dimension.
        netcdf_check(&nc_in,dim_fov = nc_in.ncid->getDim("spatial_samples_per_image").getSize());

        // FOV CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("FIELD_OF_VIEW"));

        fov_act_angles.resize(dim_fov);
        fov_ispat.resize(dim_fov*dim_detector_spec); // Contents of the file to be further processed.
        fov_ipix1.resize(dim_fov*dim_detector_spec); // Unbinned processed CKD.
        fov_ipix2.resize(dim_fov*dim_detector_spec); // Unbinned processed CKD.
        fov_weight1.resize(dim_fov*dim_detector_spec); // Unbinned processed CKD.
        fov_dims_spec.resize(dim_fov,dim_detector_spec); // Fixed content for unbinned CKD.
        dim_fov_spec_total = dim_fov*dim_detector_spec; // Folded dimension.
        fov_iel_start.resize(dim_fov);
        netcdf_check(&nc_in,grp.getVar("fov_act_angles").getVar(fov_act_angles.data()));
        fov_nfov_vp.resize(1);
        fov_nfov_vp.front() = static_cast<uint32_t>(fov_act_angles.size());
        netcdf_check(&nc_in,grp.getVar("fov_ispat").getVar(fov_ispat.data()));

        // Turn floating spatial detector coordinates into pixels and a weight factor.
        size_t iel = 0; // Folded iterator.
        for (size_t ifovpol=0 ; ifovpol<dim_fov ; ifovpol++) {
            for (size_t ispec=0 ; ispec<dim_detector_spec ; ispec++) {
                size_t ispat = static_cast<size_t>(fov_ispat[iel]);
                fov_ipix1[iel] = ispat*dim_detector_spec+ispec;
                fov_ipix2[iel] = (ispat+1)*dim_detector_spec+ispec;
                fov_weight1[iel] = ispat + 1.0 - fov_ispat[iel];
                iel++;
            }
        }

        // The fov_dims_spec are not saved in the file. CKD read from a
        // file is always unbinned, so the fov_dims_spec is always
        // dim_detector_spec for all viewports. The array is allocated
        // per FOV, because of logistic ease.
        // Fill start index array.
        for (size_t ifov=0 ; ifov<dim_fov ; ifov++) fov_iel_start[ifov] = ifov*dim_detector_spec; // This is always true for unbinned CKD.
    } // }}}
    if (lev_target > LEVEL_SWATHCAL) { // {{{

        // Swath CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("SWATH"));

        uint8_t intapply = set.swath_apply;
        swath_apply = intapply == 1;

        if (swath_apply) {
            // Swath vectors.
            swath_swathvectors.resize(dim_fov*DIM_VEC);
            netcdf_check(&nc_in,grp.getVar("swath_swathvectors").getVar(swath_swathvectors.data()));
            swath_vectorplane_normals.resize(dim_vp*DIM_VEC);
            netcdf_check(&nc_in,grp.getVar("swath_vectorplane_normals").getVar(swath_vectorplane_normals.data()));
        }
    } // }}}
    if (lev_target > LEVEL_WAVECAL) { // {{{

        // Wavelength CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("WAVELENGTH"));

        // Wavelength calibration.
        wave_spectra.resize(dim_fov*dim_detector_spec);
        wave_target.resize(dim_fov*dim_detector_spec);
        netcdf_check(&nc_in,grp.getVar("wave_spectra").getVar(wave_spectra.data()));
        // Unified spectral grid, derived from the nadirmost spectrum. It may have
        // a too high resolution. So also this one may need to be reconsidered.
        netcdf_check(&nc_in,grp.getVar("wave_target").getVar(wave_target.data()));
    } // }}}
    if (lev_target > LEVEL_RADCAL) { // {{{

        // Radiometric CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("RADIOMETRIC"));

        uint8_t intapply = set.rad_apply;
        rad_apply = intapply == 1;

        if (rad_apply) {

            // Radiometric calibration.
            rad_spectra.resize(dim_fov*dim_detector_spec);
            netcdf_check(&nc_in,grp.getVar("rad_spectra").getVar(rad_spectra.data()));
        }
    } // }}}
    return 0;

} // }}}

// TODO
// ! the check_opts fct and the write_opt fct seem not to be used anywhere.
// ! Not sure how the were envisioned to be used. Since we should not/plan not
// to write to the CKD
// ! For the moment comment out this code.
// ! At later stage the code can be deleted
// Put detector options into the CKD (to write) and warn for inconsistencies.
//int CKD::check_opts( // {{{
//    Calibration_options opt // User calibration options for current step.
//)
//{
//
//    // Execute all single-level steps.
//    if (dark_apply) handle(check_opt(opt,lev,opt_dark,LEVEL_DARKCAL,"dark"));
//    if (noise_apply) handle(check_opt(opt,lev,opt_noise,LEVEL_NOISECAL,"noise"));
//    if (nonlin_apply) handle(check_opt(opt,lev,opt_nonlin,LEVEL_NONLINCAL,"non-linearity"));
//    if (prnu_apply) handle(check_opt(opt,lev,opt_prnu,LEVEL_PRNUCAL,"PRNU"));
//    if (stray_apply) handle(check_opt(opt,lev,opt_stray,LEVEL_STRAYCAL,"straylight"));
//    handle(check_opt(opt,lev,opt_fov,LEVEL_FOVCAL,"field-of-view"));
//    if (swath_apply) handle(check_opt(opt,lev,opt_swath,LEVEL_SWATHCAL,"swath"));
//    handle(check_opt(opt,lev,opt_wave,LEVEL_WAVECAL,"wavelength"));
//    if (rad_apply) handle(check_opt(opt,lev,opt_rad,LEVEL_RADCAL,"radiometric"));
//    handle(check_opt(opt,lev,opt_pol,LEVEL_POLCAL,"polarimetric"));
//
//    return 0;
//
//} // }}}
//
//int CKD::write_opt( // {{{
//    NcGroup grp,
//    Calibration_options opt
//)
//{
//    uint8_t boolwrite; // To write any form of boolean.
//    if (lev > LEVEL_DARKCAL && dark_apply) {
//        // Forced option.
//        boolwrite = opt.dark_current?1:0;
//        netcdf_check(nc_ckd,grp.putAtt("opt_dark_current",ncUbyte,boolwrite));
//        // No optional options.
//    }
//    if (lev > LEVEL_NONLINCAL && nonlin_apply) {
//        boolwrite = opt.nonlin?1:0;
//        netcdf_check(nc_ckd,grp.putAtt("opt_nonlin",ncUbyte,boolwrite));
//        // Optional options only relevant when executing non-linearity.
//        if (opt.nonlin) {
//            netcdf_check(nc_ckd,grp.putAtt("opt_nonlin_niter",ncUint,opt.nonlin_niter));
//            netcdf_check(nc_ckd,grp.putAtt("opt_nonlin_tol",ncDouble,opt.nonlin_tol));
//        }
//    }
//    if (lev > LEVEL_STRAYCAL && stray_apply) {
//        boolwrite = opt.stray?1:0;
//        netcdf_check(nc_ckd,grp.putAtt("opt_stray",ncUbyte,boolwrite));
//        // Optional options only relevant when executing straylight.
//        if (opt.stray) {
//            netcdf_check(nc_ckd,grp.putAtt("opt_stray_van_cittert_steps",ncInt,opt.stray_van_cittert_steps));
//        }
//    }
//
//    return 0;
//
//} // }}}
//
//int CKD::check_opt( // {{{
//    Calibration_options opt_user,
//    level_t lev_user,
//    Calibration_options &opt_ref,
//    level_t lev_ref,
//    string stepname
//)
//{
//    
//    if (lev_ref < lev_user) {
//        // Warn for any inconsistencies in optional options if relevant.
//
//        // Relevant is if the reference options are high enough level that
//        // the checked options can be relevant, that both the current and
//        // the reference step execute the step and that the CKD does not
//        // skip the step.
//
//        // Non-linearity.
//        if (
//            lev_ref > LEVEL_NONLINCAL &&
//            nonlin_apply &&
//            opt_user.nonlin &&
//            opt_ref.nonlin
//        ) {
//            if (opt_user.nonlin_niter != opt_ref.nonlin_niter) writelog(log_warning,"Warning: Inconsistent non-linearity number of iterations between user optinos of current step (%d) and options when the CKD of step %s was created (%d).",opt_user.nonlin_niter,stepname.c_str(),opt_ref.nonlin_niter);
//            if (opt_user.nonlin_tol != opt_ref.nonlin_tol) writelog(log_warning,"Warning: Inconsistent non-linearity convergence criterion between user optinos of current step (%.6e) and options when the CKD of step %s was created (%.6e).",opt_user.nonlin_tol,stepname.c_str(),opt_ref.nonlin_tol);
//        }
//        // Straylight.
//        if (
//            lev_ref > LEVEL_STRAYCAL &&
//            stray_apply &&
//            opt_user.stray &&
//            opt_ref.stray
//        ) {
//            if (opt_user.stray_van_cittert_steps != opt_ref.stray_van_cittert_steps) writelog(log_warning,"Warning: Inconsistent straylight van Cittert steps between user optinos of current step (%d) and options when the CKD of step %s was created (%d).",opt_user.stray_van_cittert_steps,stepname.c_str(),opt_ref.stray_van_cittert_steps);
//        }
//    }
//    if (lev_ref == lev_user) {
//        opt_ref = opt_user; // This is the reason why opt_ref is a reference.
//    }
//    // if lev_ref is greater than lev_user, nothing needs to be done.
//
//    // If lev_user is L1B, there are only checks, no writes.
//
//    return 0;

//} // }}}

} // namespace tango
