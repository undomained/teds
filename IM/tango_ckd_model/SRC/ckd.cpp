#include "header.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"

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

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_DARKCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("DARK"));

        // Define diagnostic variable.
        vector<double> dark_chi2;

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("dark_skip").getVar(&intskip));

        dark_skip = intskip == 1;

        if (!dark_skip) {

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

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_NOISECAL) { // {{{

        // Noise CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("NOISE"));

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("noise_skip").getVar(&intskip));

        noise_skip = intskip == 1;

        if (!noise_skip) {

            // Shape the noise CKD.
            noise_g.resize(npix);
            noise_n.resize(npix);

            // Read the noise CKD.
            netcdf_check(&nc_in,grp.getVar("noise_g").getVar(noise_g.data()));
            netcdf_check(&nc_in,grp.getVar("noise_n").getVar(noise_n.data()));
        }

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_NONLINCAL) { // {{{

        // Non-linearity CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("NON_LINEARITY"));

        // Define diagnostic variables.
        vector<double> nonlin_lin_slope;
        vector<double> nonlin_lin_chi2;
        vector<double> nonlin_chi2;

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("nonlin_skip").getVar(&intskip));

        nonlin_skip = intskip == 1;

        if (!nonlin_skip) {

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

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_PRNUCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("PRNU"));

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("prnu_skip").getVar(&intskip));

        prnu_skip = intskip == 1;

        if (!prnu_skip) {
            // Pixel response non-uniformity.
            prnu_prnu.resize(npix);
            netcdf_check(&nc_in,grp.getVar("prnu_prnu").getVar(prnu_prnu.data()));
        }

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_STRAYCAL) { // {{{

        // Dark CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("STRAYLIGHT"));

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("stray_skip").getVar(&intskip));

        stray_skip = intskip == 1;

        if (!stray_skip) {
            netcdf_check(&nc_in,stray_kernel_n_rows = grp.getDim("stray_kernel_n_rows").getSize());
            netcdf_check(&nc_in,stray_kernel_n_cols = grp.getDim("stray_kernel_n_cols").getSize());
            // Need to divide by 2 because of complex numbers
            netcdf_check(&nc_in,stray_kernel_fft_size = grp.getDim("stray_kernel_fft_size").getSize()/2);
            stray_kernel_fft.resize(stray_kernel_fft_size);
            netcdf_check(&nc_in,grp.getVar("stray_kernel_fft").getVar(stray_kernel_fft.data()));
            stray_moving_kernel_fft.resize(stray_kernel_fft_size);
            netcdf_check(&nc_in,grp.getVar("stray_moving_kernel_fft").getVar(stray_moving_kernel_fft.data()));
            netcdf_check(&nc_in,stray_transformed_n_rows = grp.getDim("stray_transformed_n_rows").getSize());
            netcdf_check(&nc_in,stray_transformed_n_cols = grp.getDim("stray_transformed_n_cols").getSize());
            stray_transform_indices.resize(stray_transformed_n_rows*stray_transformed_n_cols);
            netcdf_check(&nc_in,grp.getVar("stray_transform_indices").getVar(stray_transform_indices.data()));
            stray_transform_deltas.resize(4*stray_transformed_n_rows*stray_transformed_n_cols);
            netcdf_check(&nc_in,grp.getVar("stray_transform_deltas").getVar(stray_transform_deltas.data()));
            netcdf_check(&nc_in,grp.getVar("stray_eta").getVar(&stray_eta));

            if (const auto g { grp.getGroup("interpolating") }; !g.isNull()) {
                stray.n_kernels =
                  static_cast<int>(g.getDim("kernels").getSize());
                stray.n_spatial =
                  static_cast<int>(g.getDim("spatial").getSize());
                stray.n_spectral =
                  static_cast<int>(g.getDim("spectral").getSize());
                stray.kernel_rows.resize(stray.n_kernels);
                g.getVar("kernel_rows").getVar(stray.kernel_rows.data());
                stray.kernel_cols.resize(stray.n_kernels);
                g.getVar("kernel_cols").getVar(stray.kernel_cols.data());
                stray.kernel_fft_sizes.resize(stray.n_kernels);
                g.getVar("kernel_fft_sizes")
                  .getVar(stray.kernel_fft_sizes.data());
                stray.kernels_fft.resize(stray.n_kernels);
                for (int i {}; i < stray.n_kernels; ++i) {
                    stray.kernels_fft[i].resize(stray.kernel_fft_sizes[i]);
                }
                std::vector<size_t> start { 0 };
                for (int i_kernel {}; i_kernel < stray.n_kernels; ++i_kernel) {
                    // Need to multiply by 2 because of complex numbers
                    std::vector<size_t> count {
                        2 * static_cast<size_t>(
                          stray.kernel_fft_sizes.at(i_kernel))
                    };
                    g.getVar("kernels_fft")
                      .getVar(start,
                              count,
                              stray.kernels_fft.at(i_kernel).data());
                    start.front() += count.front();
                }
                stray.eta.resize(stray.n_spatial * stray.n_spectral);
                g.getVar("eta").getVar(stray.eta.data());
                stray.weights.resize(
                  stray.n_kernels,
                  std::vector<double>(stray.n_spatial * stray.n_spectral));
                std::vector<double> buf(stray.n_kernels * stray.n_spatial
                                        * stray.n_spectral);
                g.getVar("weights").getVar(buf.data());
                for (int i {}; i < stray.n_kernels; ++i) {
                    for (int j {}; j < stray.n_spatial * stray.n_spectral;
                         ++j) {
                        stray.weights[i][j] =
                          buf[i * stray.n_spatial * stray.n_spectral + j];
                    }
                }
                stray.edges.resize(4 * stray.n_kernels);
                g.getVar("edges").getVar(stray.edges.data());
            }
        }

        // Copy this input to the output.
        handle(writestep());

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

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_SWATHCAL) { // {{{

        // Swath CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("SWATH"));

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("swath_skip").getVar(&intskip));

        swath_skip = intskip == 1;

        if (!swath_skip) {
            // Swath vectors.
            swath_swathvectors.resize(dim_fov*DIM_VEC);
            netcdf_check(&nc_in,grp.getVar("swath_swathvectors").getVar(swath_swathvectors.data()));
            swath_vectorplane_normals.resize(dim_vp*DIM_VEC);
            netcdf_check(&nc_in,grp.getVar("swath_vectorplane_normals").getVar(swath_vectorplane_normals.data()));
        }

        // Copy this input to the output.
        handle(writestep());

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

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    if (lev_target > LEVEL_RADCAL) { // {{{

        // Radiometric CKD group.
        NcGroup grp;
        netcdf_check(&nc_in,grp = nc_in.ncid->getGroup("RADIOMETRIC"));

        uint8_t intskip;
        netcdf_check(&nc_in,grp.getVar("rad_skip").getVar(&intskip));

        rad_skip = intskip == 1;

        if (!rad_skip) {

            // Radiometric calibration.
            rad_spectra.resize(dim_fov*dim_detector_spec);
            netcdf_check(&nc_in,grp.getVar("rad_spectra").getVar(rad_spectra.data()));
        }

        // Copy this input to the output.
        handle(writestep());

    } // }}}
    return 0;

} // }}}

// Writes the current step into the CKD making the total CKD one step more
// mature (increasing its level by one).
int CKD::writestep( // {{{
)
{

    // This routine only does something if the level is a CKD-generation step.
    // It writes all the relevant CKD including the diagnostic CKD. Pointers to
    // diagnostic CKD are in the CKD structure and during this routine their
    // targets should be alive, but later on, the targets are probably gone.
    if (write_ckd) {
        vector<uint8_t> mask_values = {0,1}; // For skipping flags or masks.
        NcVar var; // Variable that is directly filled.
        if (lev == LEVEL_DIMCAL) { // {{{

            // Basic dimensions.
            netcdf_check(nc_ckd,dimid_pol = nc_ckd->ncid->addDim("polarization_directions",DIM_POL));
            netcdf_check(nc_ckd,dimid_vec = nc_ckd->ncid->addDim("vector_elements",DIM_VEC));
            netcdf_check(nc_ckd,dimid_detector_spec = nc_ckd->ncid->addDim("spectral_detector_pixels",dim_detector_spec));
            netcdf_check(nc_ckd,dimid_detector_spat = nc_ckd->ncid->addDim("spatial_detector_pixels",dim_detector_spat));
            netcdf_check(nc_ckd,dimid_vp = nc_ckd->ncid->addDim("number_of_views",dim_vp));

            // Define the pixel mask.
            vector<NcDim> dims_mask = {dimid_detector_spat,dimid_detector_spec};
            netcdf_check(nc_ckd,var_mask = nc_ckd->ncid->addVar("mask",ncUbyte,dims_mask));
            netcdf_check(nc_ckd,var_mask.putAtt("long_name","detector pixel mask"));
            netcdf_check(nc_ckd,var_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var_mask.putAtt("flag_meanings","Good Bad"));

            // Define the viewport mask.
            netcdf_check(nc_ckd,var_vp_mask = nc_ckd->ncid->addVar("vp_mask",ncUbyte,dimid_vp));
            netcdf_check(nc_ckd,var_vp_mask.putAtt("long_name","viewport mask"));
            netcdf_check(nc_ckd,var_vp_mask.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var_vp_mask.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var_vp_mask.putAtt("comment","Exclusion only used in test environment"));

        } // }}}
        if (lev == LEVEL_DARKCAL) { // {{{

            // Dark CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("DARK"));

            // Dark skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("dark_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","dark calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = dark_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!dark_skip) {

                // Create order dimension for temperature.
                NcDim dimid_dark_order;
                netcdf_check(nc_ckd,dimid_dark_order = grp.addDim("dark_number_of_coefficients",dim_dark_order));

                // Create variables for dark.
                // Dark offset.
                netcdf_check(nc_ckd,var = grp.addVar("dark_offset",ncDouble,{dimid_detector_spat,dimid_detector_spec,dimid_dark_order}));
                netcdf_check(nc_ckd,var.putAtt("long_name","b-spline coefficients for signal offset as function of temperature"));
                vector<double> knots = {dark_nominal_temperature,dark_nominal_temperature+1.0};
                netcdf_check(nc_ckd,var.putAtt("knots",ncDouble,knots.size(),knots.data()));
                netcdf_check(nc_ckd,var.putAtt("knot_units","K"));
                netcdf_check(nc_ckd,var.putAtt("comment","basis set continues as polynomials outside knot domain"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putVar(dark_offset.data()));

                // Dark current.
                netcdf_check(nc_ckd,var = grp.addVar("dark_current",ncDouble,{dimid_detector_spat,dimid_detector_spec,dimid_dark_order}));
                netcdf_check(nc_ckd,var.putAtt("long_name","b-spline coefficients for dark current as function of temperature"));
                // Recycle knots vector from dark offset.
                netcdf_check(nc_ckd,var.putAtt("knots",ncDouble,knots.size(),knots.data()));
                netcdf_check(nc_ckd,var.putAtt("knot_units","K"));
                netcdf_check(nc_ckd,var.putAtt("comment","basis set continues as polynomials outside knot domain"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts / second"));
                netcdf_check(nc_ckd,var.putVar(dark_current.data()));

                // Dark retrieval residual chi squared.
                netcdf_check(nc_ckd,var = grp.addVar("dark_chi2",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","residual chi squared of dark offset and current fit"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putAtt("comment","not scaled with noise"));
                netcdf_check(nc_ckd,var.putVar(diag_dark_chi2));

                // Dark nominal temperature.
                netcdf_check(nc_ckd,var = grp.addVar("dark_nominal_temperature",ncDouble));
                netcdf_check(nc_ckd,var.putAtt("long_name","detector temperature at nominal conditions"));
                netcdf_check(nc_ckd,var.putAtt("units","K"));
                netcdf_check(nc_ckd,var.putVar(&dark_nominal_temperature));

                handle(write_opt(grp,opt_dark));

            }

        } // }}}
        if (lev == LEVEL_NOISECAL) { // {{{

            // Noise CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("NOISE"));

            // Noise skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("noise_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","noise calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = noise_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!noise_skip) {

                // Write noise CKD.
                // Noise g-paramter.
                netcdf_check(nc_ckd,var = grp.addVar("noise_g",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","parameter g in noise model"));
                netcdf_check(nc_ckd,var.putAtt("noise_model","noise = sqrt(g*signal + n^2)"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putAtt("related_parameter","noise_n"));
                netcdf_check(nc_ckd,var.putVar(noise_g.data()));

                // Noise n-paramter.
                netcdf_check(nc_ckd,var = grp.addVar("noise_n",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","parameter n in noise model"));
                netcdf_check(nc_ckd,var.putAtt("noise_model","noise = sqrt(g*signal + n^2)"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putAtt("related_parameter","noise_g"));
                netcdf_check(nc_ckd,var.putVar(noise_n.data()));

                handle(write_opt(grp,opt_noise));

            }

        } // }}}
        if (lev == LEVEL_NONLINCAL) { // {{{

            // Non-linearity CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("NON_LINEARITY"));

            // Non-linearity skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("nonlin_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","non-linearity calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = nonlin_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!nonlin_skip) {
                // These dimensions identifiers are only needed to write this output.
                // The series dimension is only for diagnostic CKD, so its dimension size is
                // even a non-owned pointer.
                NcDim dimid_nonlin_exptime;
                netcdf_check(nc_ckd,dimid_nonlin_exptime = grp.addDim("nonlin_number_of_exposure_times",dim_nonlin_exptime));
                NcDim dimid_nonlin_spline;
                netcdf_check(nc_ckd,dimid_nonlin_spline = grp.addDim("nonlin_number_of_coefficients",dim_nonlin_spline));
                NcDim dimid_nonlin_knot;
                netcdf_check(nc_ckd,dimid_nonlin_knot = grp.addDim("nonlin_number_of_knots",dim_nonlin_knot));

                // Linear fit slope.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_lin_slope",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","fitted slope in linear regime"));
                netcdf_check(nc_ckd,var.putAtt("units","0.001 detector counts / electron"));
                netcdf_check(nc_ckd,var.putVar(diag_nonlin_lin_slope));

                // Linear fit residual chi square.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_lin_chi2",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","residual chi squared of linear fit through linear regime"));
                netcdf_check(nc_ckd,var.putAtt("units","1"));
                netcdf_check(nc_ckd,var.putVar(diag_nonlin_lin_chi2));

                // Order (scalar).
                int intorder = (int) nonlin_order;
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_order",ncInt));
                netcdf_check(nc_ckd,var.putAtt("long_name","b-spline order of non-linearity fit"));
                netcdf_check(nc_ckd,var.putAtt("definition","parameter n, so polynomial order is one lower"));
                netcdf_check(nc_ckd,var.putVar(&intorder));

                // Knots.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_knots",ncDouble,dimid_nonlin_knot));
                netcdf_check(nc_ckd,var.putAtt("long_name","scaled b-spline knots of non-linearity fit"));
                netcdf_check(nc_ckd,var.putAtt("units","1"));
                netcdf_check(nc_ckd,var.putAtt("related_parameter","nonlin_signal_scale"));
                netcdf_check(nc_ckd,var.putVar(nonlin_knots.data()));

                // Exposure times.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_exptimes",ncDouble,dimid_nonlin_exptime));
                netcdf_check(nc_ckd,var.putAtt("long_name","exposure times for which non-linearity correction is derived"));
                netcdf_check(nc_ckd,var.putAtt("units","s"));
                netcdf_check(nc_ckd,var.putVar(nonlin_exptimes.data()));

                // Signal scales.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_signal_scale",ncDouble,{dimid_detector_spat,dimid_detector_spec,dimid_nonlin_exptime}));
                netcdf_check(nc_ckd,var.putAtt("long_name","scaling factor of the b-spline knots"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putAtt("related_parameter","nonlin_knots"));
                netcdf_check(nc_ckd,var.putVar(nonlin_signal_scale.data()));

                // Fit results.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_fit",ncDouble,{dimid_detector_spat,dimid_detector_spec,dimid_nonlin_exptime,dimid_nonlin_spline}));
                netcdf_check(nc_ckd,var.putAtt("long_name","b-spline coefficients for non-linearity curve"));
                netcdf_check(nc_ckd,var.putAtt("units","detector counts"));
                netcdf_check(nc_ckd,var.putVar(nonlin_fit.data()));

                // Fit residual chi square.
                netcdf_check(nc_ckd,var = grp.addVar("nonlin_chi2",ncDouble,{dimid_detector_spat,dimid_detector_spec,dimid_nonlin_exptime}));
                netcdf_check(nc_ckd,var.putAtt("long_name","residual chi squared of non-linearity curve fit"));
                netcdf_check(nc_ckd,var.putAtt("units","1"));
                netcdf_check(nc_ckd,var.putVar(diag_nonlin_chi2));

                handle(write_opt(grp,opt_nonlin));

            }

        } // }}}
        if (lev == LEVEL_PRNUCAL) { // {{{

            // PRNU CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("PRNU"));

            // PRNU skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("prnu_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","PRNU calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = prnu_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!prnu_skip) {
                // Pixel-response non-uniformity.
                netcdf_check(nc_ckd,var = grp.addVar("prnu_prnu",ncDouble,{dimid_detector_spat,dimid_detector_spec}));
                netcdf_check(nc_ckd,var.putAtt("long_name","relative photon sensitivity"));
                netcdf_check(nc_ckd,var.putAtt("units","1"));
                netcdf_check(nc_ckd,var.putVar(prnu_prnu.data()));

                handle(write_opt(grp,opt_prnu));
            }

        } // }}}
        if (lev == LEVEL_STRAYCAL) { // {{{

            // Straylight CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("STRAYLIGHT"));

            // Straylight skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("stray_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","Straylight correction execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = stray_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!stray_skip) {
                NcDim stray_kernel_n_rows_id;
                netcdf_check(nc_ckd, stray_kernel_n_rows_id = grp.addDim("stray_kernel_n_rows", stray_kernel_n_rows));
                NcDim stray_kernel_n_cols_id;
                netcdf_check(nc_ckd, stray_kernel_n_cols_id = grp.addDim("stray_kernel_n_cols", stray_kernel_n_cols));
                NcDim stray_kernel_fft_size_id;
                // Need to multiply by 2 because of complex numbers
                netcdf_check(nc_ckd, stray_kernel_fft_size_id = grp.addDim("stray_kernel_fft_size", 2*stray_kernel_fft_size));
                netcdf_check(nc_ckd, grp.addVar("stray_kernel_fft", ncDouble, stray_kernel_fft_size_id).putVar(stray_kernel_fft.data()));
                netcdf_check(nc_ckd, grp.addVar("stray_moving_kernel_fft", ncDouble, stray_kernel_fft_size_id).putVar(stray_moving_kernel_fft.data()));
                NcDim stray_transformed_n_rows_id;
                netcdf_check(nc_ckd, stray_transformed_n_rows_id = grp.addDim("stray_transformed_n_rows", stray_transformed_n_rows));
                NcDim stray_transformed_n_cols_id;
                netcdf_check(nc_ckd, stray_transformed_n_cols_id = grp.addDim("stray_transformed_n_cols", stray_transformed_n_cols));
                vector<NcDim> stray_transformed_dims { stray_transformed_n_rows_id, stray_transformed_n_cols_id };
                netcdf_check(nc_ckd, grp.addVar("stray_transform_indices", ncInt, stray_transformed_dims).putVar(stray_transform_indices.data()));
                NcDim stray_transform_deltas_size_id;
                netcdf_check(nc_ckd, stray_transform_deltas_size_id = grp.addDim("stray_transform_deltas_size", 4*stray_transformed_n_rows*stray_transformed_n_cols));
                netcdf_check(nc_ckd, grp.addVar("stray_transform_deltas", ncDouble, stray_transform_deltas_size_id).putVar(stray_transform_deltas.data()));
                netcdf_check(nc_ckd, grp.addVar("stray_eta", ncDouble).putVar(&stray_eta));
                const uint8_t dry_run { static_cast<uint8_t>(stray_dry_run) };
                netcdf_check(nc_ckd, grp.addVar("stray_dry_run", ncUbyte).putVar(&dry_run));

                handle(write_opt(grp,opt_stray));
            }

        } // }}}
        if (lev == LEVEL_FOVCAL) { // {{{

            // New dimension: fov.
            // This dimensions must be defined in the root group, because
            // it is used by different steps.
            netcdf_check(nc_ckd,dimid_fov = nc_ckd->ncid->addDim("spatial_samples_per_image",dim_fov));

            // FOV CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("FIELD_OF_VIEW"));

            // Create variables for FOV.
            // Number of spatial samples per viewport.
            netcdf_check(nc_ckd,var = grp.addVar("fov_nfov_vp",ncUint,dimid_vp));
            netcdf_check(nc_ckd,var.putAtt("long_name","number of spatial samples per viewport"));
            netcdf_check(nc_ckd,var.putVar(fov_nfov_vp.data()));

            // Across-track angles per viewport.
            netcdf_check(nc_ckd,var = grp.addVar("fov_act_angles",ncDouble,dimid_fov));
            netcdf_check(nc_ckd,var.putAtt("long_name","across-track rotation stage position per spatial sample"));
            netcdf_check(nc_ckd,var.putAtt("units","radians"));
            netcdf_check(nc_ckd,var.putVar(fov_act_angles.data()));

            vector<NcDim> ispatdims = {dimid_fov,dimid_pol,dimid_detector_spec};

            // Positions on the detector of the spectra.
            netcdf_check(nc_ckd,var = grp.addVar("fov_ispat",ncDouble,ispatdims));
            netcdf_check(nc_ckd,var.putAtt("long_name","floating point spatial detector pixel indices of the spectra"));
            netcdf_check(nc_ckd,var.putAtt("units","unbinned detector pixels"));
            netcdf_check(nc_ckd,var.putVar(fov_ispat.data()));

            handle(write_opt(grp,opt_fov));

        } // }}}
        if (lev == LEVEL_SWATHCAL) { // {{{

            // Swath CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("SWATH"));

            // Swath skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("swath_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","swath vector calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = swath_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!swath_skip) {
                // Swath vectors.
                vector<NcDim> dims_swathvectors = {dimid_fov,dimid_vec};
                netcdf_check(nc_ckd,var = grp.addVar("swath_swathvectors",ncDouble,dims_swathvectors));
                netcdf_check(nc_ckd,var.putAtt("long_name","normalized pointing vectors in spacecraft coordinates"));
                netcdf_check(nc_ckd,var.putAtt("axes","x = flying direction, y = to the right, z = down"));
                netcdf_check(nc_ckd,var.putVar(swath_swathvectors.data()));
                vector<NcDim> dims_vectorplane_normals = {dimid_vp,dimid_vec};
                netcdf_check(nc_ckd,var = grp.addVar("swath_vectorplane_normals",ncDouble,dims_vectorplane_normals));
                netcdf_check(nc_ckd,var.putAtt("long_name","Normal vectors to plane through swath vectors of viewport."));
                netcdf_check(nc_ckd,var.putAtt("axes","x = flying direction, y = to the right, z = down"));
                netcdf_check(nc_ckd,var.putVar(swath_vectorplane_normals.data()));

                handle(write_opt(grp,opt_swath));
            }

        } // }}}
        if (lev == LEVEL_WAVECAL) { // {{{

            // Wavelength CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("WAVELENGTH"));

            // Wavelength calibration. This is just a map of wavelengths for the extracted
            // spectra.

            // Wavelengths per spectrum.
            vector<NcDim> dims_wave_spectra = {dimid_fov,dimid_pol,dimid_detector_spec};
            netcdf_check(nc_ckd,var = grp.addVar("wave_spectra",ncDouble,dims_wave_spectra));
            netcdf_check(nc_ckd,var.putAtt("long_name","wavelengths of S and P spectra"));
            netcdf_check(nc_ckd,var.putAtt("units","nm"));
            netcdf_check(nc_ckd,var.putVar(wave_spectra.data()));

            // Interpolation wavelengths for demodulation.
            vector<NcDim> dims_wave_target = {dimid_fov,dimid_detector_spec};
            netcdf_check(nc_ckd,var = grp.addVar("wave_target",ncDouble,dims_wave_target));
            netcdf_check(nc_ckd,var.putAtt("long_name","wavelengths of modulation spectra"));
            netcdf_check(nc_ckd,var.putAtt("units","nm"));
            netcdf_check(nc_ckd,var.putVar(wave_target.data()));

            handle(write_opt(grp,opt_wave));

        } // }}}
        if (lev == LEVEL_RADCAL) { // {{{

            // Radiometric CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("RADIOMETRIC"));

            // Radiometric calibration skip flag.
            netcdf_check(nc_ckd,var = grp.addVar("rad_skip",ncUbyte));
            netcdf_check(nc_ckd,var.putAtt("long_name","radiometric calibration calibration execution flag"));
            netcdf_check(nc_ckd,var.putAtt("flag_values",ncUbyte,mask_values.size(),mask_values.data()));
            netcdf_check(nc_ckd,var.putAtt("flag_meanings","Included Excluded"));
            netcdf_check(nc_ckd,var.putAtt("comment","Exclusion only used in test environment"));
            uint8_t intskip = rad_skip?1:0;
            netcdf_check(nc_ckd,var.putVar(&intskip));

            if (!rad_skip) {
                // Radiometric calibration data on spectrum perspective.
                vector<NcDim> dims_rad_spectra = {dimid_fov,dimid_pol,dimid_detector_spec};
                netcdf_check(nc_ckd,var = grp.addVar("rad_spectra",ncDouble,dims_rad_spectra));
                netcdf_check(nc_ckd,var.putAtt("long_name","radiometric calibration constant"));
                netcdf_check(nc_ckd,var.putAtt("units","J.m-2.sr-1.um-1 / detector count"));
                netcdf_check(nc_ckd,var.putVar(rad_spectra.data()));

                handle(write_opt(grp,opt_rad));
            }

        } // }}}
        if (lev == LEVEL_POLCAL) { // {{{

            // Dark CKD group.
            NcGroup grp;
            netcdf_check(nc_ckd,grp = nc_ckd->ncid->addGroup("POLARIMETRIC"));

            // Polarimetric calibration sampled on spectral pixels.
            vector<NcDim> dims_pol = {dimid_fov,dimid_detector_spec};

            // Mueller matrix ratio for Q.
            netcdf_check(nc_ckd,var = grp.addVar("pol_m_q",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","modulation sensitivy to q-polarization"));
            netcdf_check(nc_ckd,var.putAtt("units","1"));
            netcdf_check(nc_ckd,var.putVar(pol_m_q.data()));

            // Mueller matrix ratio for U.
            netcdf_check(nc_ckd,var = grp.addVar("pol_m_u",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","modulation sensitivy to u-polarization"));
            netcdf_check(nc_ckd,var.putAtt("units","1"));
            netcdf_check(nc_ckd,var.putVar(pol_m_u.data()));

            // Mueller element for telescope polarization.
            netcdf_check(nc_ckd,var = grp.addVar("pol_m_t",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","telescope polarization"));
            netcdf_check(nc_ckd,var.putAtt("units","1"));
            netcdf_check(nc_ckd,var.putVar(pol_m_t.data()));

            // Retardance.
            netcdf_check(nc_ckd,var = grp.addVar("pol_delta",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","retardance"));
            netcdf_check(nc_ckd,var.putAtt("units","nm"));
            netcdf_check(nc_ckd,var.putVar(diag_pol_delta));

            // Modulation efficiency axis A.
            netcdf_check(nc_ckd,var = grp.addVar("pol_eff_a",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","modulation efficiency on first axis"));
            netcdf_check(nc_ckd,var.putAtt("units","1"));
            netcdf_check(nc_ckd,var.putAtt("axis","q-polarization for zero tilt"));
            netcdf_check(nc_ckd,var.putVar(diag_pol_eff_a));

            // Modulation efficiency axis B.
            netcdf_check(nc_ckd,var = grp.addVar("pol_eff_b",ncDouble,dims_pol));
            netcdf_check(nc_ckd,var.putAtt("long_name","modulation efficiency on second axis"));
            netcdf_check(nc_ckd,var.putAtt("units","1"));
            netcdf_check(nc_ckd,var.putAtt("axis","u-polarization for zero tilt"));
            netcdf_check(nc_ckd,var.putVar(diag_pol_eff_b));

            // Tilt angle.
            netcdf_check(nc_ckd,var = grp.addVar("pol_tilt",ncDouble,dimid_fov));
            netcdf_check(nc_ckd,var.putAtt("long_name","tilt between axis system a/b and q/u"));
            netcdf_check(nc_ckd,var.putAtt("units","degrees"));
            netcdf_check(nc_ckd,var.putVar(diag_pol_tilt));

            handle(write_opt(grp,opt_pol));

        } // }}}

        // Update pixel mask, by any of the processor steps that do something with it.
        if (lev == LEVEL_DIMCAL || lev == LEVEL_NONLINCAL || lev == LEVEL_DARKCAL || lev == LEVEL_NOISECAL || lev == LEVEL_PRNUCAL) {

            // Convert boolean to unsigned 8-bit integers.
            vector<uint8_t> writemask(npix);
            for (size_t ipix=0 ; ipix<npix ; ipix++) writemask[ipix] = mask[ipix]?1:0;
            netcdf_check(nc_ckd,var_mask.putVar(writemask.data()));

        }

        // Update viewport mask. This can be done after any step, though
        // it is logical only to do so at a viewport-related step.
        // The point is that the viewport mask is updated by the parent
        // class.
        vector<uint8_t> vp_writemask(dim_vp);
        for (size_t ivp=0 ; ivp<dim_vp ; ivp++) vp_writemask[ivp] = vp_mask[ivp]?1:0;
        netcdf_check(nc_ckd,var_vp_mask.putVar(vp_writemask.data()));

        // Synchronize, so that an error in a later step does not destroy
        // all output.
        netcdf_check(nc_ckd,nc_ckd->ncid->sync());

    }

    // Progress one level.
    lev = (level_t) ((int)lev + 1);

    return 0;

} // }}}

// Put detector options into the CKD (to write) and warn for inconsistencies.
int CKD::check_opts( // {{{
    Calibration_options opt // User calibration options for current step.
)
{

    // Execute all single-level steps.
    if (!dark_skip) handle(check_opt(opt,lev,opt_dark,LEVEL_DARKCAL,"dark"));
    if (!noise_skip) handle(check_opt(opt,lev,opt_noise,LEVEL_NOISECAL,"noise"));
    if (!nonlin_skip) handle(check_opt(opt,lev,opt_nonlin,LEVEL_NONLINCAL,"non-linearity"));
    if (!prnu_skip) handle(check_opt(opt,lev,opt_prnu,LEVEL_PRNUCAL,"PRNU"));
    if (!stray_skip) handle(check_opt(opt,lev,opt_stray,LEVEL_STRAYCAL,"straylight"));
    handle(check_opt(opt,lev,opt_fov,LEVEL_FOVCAL,"field-of-view"));
    if (!swath_skip) handle(check_opt(opt,lev,opt_swath,LEVEL_SWATHCAL,"swath"));
    handle(check_opt(opt,lev,opt_wave,LEVEL_WAVECAL,"wavelength"));
    if (!rad_skip) handle(check_opt(opt,lev,opt_rad,LEVEL_RADCAL,"radiometric"));
    handle(check_opt(opt,lev,opt_pol,LEVEL_POLCAL,"polarimetric"));

    return 0;

} // }}}

int CKD::write_opt( // {{{
    NcGroup grp,
    Calibration_options opt
)
{
    uint8_t boolwrite; // To write any form of boolean.
    if (lev > LEVEL_DARKCAL && !dark_skip) {
        // Forced option.
        boolwrite = opt.dark_current?1:0;
        netcdf_check(nc_ckd,grp.putAtt("opt_dark_current",ncUbyte,boolwrite));
        // No optional options.
    }
    if (lev > LEVEL_NONLINCAL && !nonlin_skip) {
        boolwrite = opt.nonlin?1:0;
        netcdf_check(nc_ckd,grp.putAtt("opt_nonlin",ncUbyte,boolwrite));
        // Optional options only relevant when executing non-linearity.
        if (opt.nonlin) {
            netcdf_check(nc_ckd,grp.putAtt("opt_nonlin_niter",ncUint,opt.nonlin_niter));
            netcdf_check(nc_ckd,grp.putAtt("opt_nonlin_tol",ncDouble,opt.nonlin_tol));
        }
    }
    if (lev > LEVEL_STRAYCAL && !stray_skip) {
        boolwrite = opt.stray?1:0;
        netcdf_check(nc_ckd,grp.putAtt("opt_stray",ncUbyte,boolwrite));
        // Optional options only relevant when executing straylight.
        if (opt.stray) {
            netcdf_check(nc_ckd,grp.putAtt("opt_stray_van_cittert_steps",ncInt,opt.stray_van_cittert_steps));
        }
    }

    return 0;

} // }}}

int CKD::check_opt( // {{{
    Calibration_options opt_user,
    level_t lev_user,
    Calibration_options &opt_ref,
    level_t lev_ref,
    string stepname
)
{
    
    if (lev_ref < lev_user) {
        // Warn for any inconsistencies in optional options if relevant.

        // Relevant is if the reference options are high enough level that
        // the checked options can be relevant, that both the current and
        // the reference step execute the step and that the CKD does not
        // skip the step.

        // Non-linearity.
        if (
            lev_ref > LEVEL_NONLINCAL &&
            !nonlin_skip &&
            opt_user.nonlin &&
            opt_ref.nonlin
        ) {
            if (opt_user.nonlin_niter != opt_ref.nonlin_niter) writelog(log_warning,"Warning: Inconsistent non-linearity number of iterations between user optinos of current step (%d) and options when the CKD of step %s was created (%d).",opt_user.nonlin_niter,stepname.c_str(),opt_ref.nonlin_niter);
            if (opt_user.nonlin_tol != opt_ref.nonlin_tol) writelog(log_warning,"Warning: Inconsistent non-linearity convergence criterion between user optinos of current step (%.6e) and options when the CKD of step %s was created (%.6e).",opt_user.nonlin_tol,stepname.c_str(),opt_ref.nonlin_tol);
        }
        // Straylight.
        if (
            lev_ref > LEVEL_STRAYCAL &&
            !stray_skip &&
            opt_user.stray &&
            opt_ref.stray
        ) {
            if (opt_user.stray_van_cittert_steps != opt_ref.stray_van_cittert_steps) writelog(log_warning,"Warning: Inconsistent straylight van Cittert steps between user optinos of current step (%d) and options when the CKD of step %s was created (%d).",opt_user.stray_van_cittert_steps,stepname.c_str(),opt_ref.stray_van_cittert_steps);
        }
    }
    if (lev_ref == lev_user) {
        opt_ref = opt_user; // This is the reason why opt_ref is a reference.
    }
    // if lev_ref is greater than lev_user, nothing needs to be done.

    // If lev_user is L1B, there are only checks, no writes.

    return 0;

} // }}}

