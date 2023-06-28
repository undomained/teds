// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"
#include "ckd.h"
#include "binningtable.h"

// Binning table file. All binning tables are included in one file.
// That file governs the pre-binning, which is applied for all non-trivial
// binning tables. This means that all binning tables used at the same time
// use the same prebinning, except for possibly the trivial binning table
// that is just 'no binning at all, even no prebinning'. A binning table
// with just pre-binning should be defined with a table that does no
// furtherm binning.

// Constructor.
Binningtable_file::Binningtable_file( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Binningtable_file::~Binningtable_file() {} // Destructor.

int Binningtable_file::read( // {{{
    string &binningtable_filename, // File with all the used binning tables.
    CKD *ckd // Pointer to the general CKD.
)
{

    // Open the NetCDF file. Let us hope we can keep it open. It will be
    // closed when the L1A manager dies.
    nc = make_unique<NetCDF_object>(this);
    handle(nc->open(binningtable_filename,NcFile::read));
    // Read dimensions in main group.
    netcdf_check(nc,dim_binningtable_spat = nc->ncid->getDim("row").getSize());
    netcdf_check(nc,dim_binningtable_spec = nc->ncid->getDim("column").getSize());
    // Verify that the pre-binning is with integer factor.
    prebinning_spat = ckd->dim_detector_spat / dim_binningtable_spat; // This must be an exact division.
    check_error(prebinning_spat * dim_binningtable_spat != ckd->dim_detector_spat,"Error: Prebinning factor in spatial (row) dimension is not integral: %zu to %zu.",ckd->dim_detector_spat,dim_binningtable_spat);
    prebinning_spec = ckd->dim_detector_spec / dim_binningtable_spec; // This must be an exact division.
    check_error(prebinning_spec * dim_binningtable_spec != ckd->dim_detector_spec,"Error: Prebinning factor in spectral (column) dimension is not integral: %zu to %zu.",ckd->dim_detector_spec,dim_binningtable_spec);
    // Leave the file open, it can be used by the binning tables themselves.

    return 0;

} // }}}

// Constructor.
Binningtable::Binningtable( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Binningtable::~Binningtable() {} // Destructor.

// Getter.
size_t Binningtable::getBinsize( // {{{
    size_t ipix
)
{
    // Return size of the bin. If the binning table is trivial. Everything
    // has size one.
    if (trivial) return 1;
    else return binsizes[ipix];
} // }}}

// Binning table read routine.
int Binningtable::read( // {{{
    Binningtable_file *bin_file, // Binning table file.
    uint8_t id, // Mask identifier.
    CKD *ckd // Pointer to the general CKD.
)
{

    check_error(ckd->lev <= LEVEL_DIMCAL,"Error: No L1A image can be managed since the CKD does not yet contain the detector shape.");
    npix_unbinned = ckd->npix; // So, the CKD need not be dragged to the private routine average_ckd.

    if (id == 0) {
        // Trivial table.
        trivial = true;
        binned_ckd = ckd;
        npix = npix_unbinned;
    } else {
        trivial = false;
        NetCDF_object *nc = bin_file->nc.get(); // Short-access.
        NcGroup grp;
        // string groupname = format("Table_%2.2i",id);
        string groupname = format("Table_{:02d}",id);
        netcdf_check(nc,grp = nc->ncid->getGroup(groupname.c_str()));
        vector<uint32_t> binning_table_read(bin_file->dim_binningtable_spat*bin_file->dim_binningtable_spec);
        netcdf_check(nc,grp.getVar("binning_table").getVar(binning_table_read.data())); // From file, that is after pre-binning.
        pixelpointer.resize(npix_unbinned); // This is the pixel pointer of the combined binning.
        // Count the sizes of each bin (used a lot during averaging).
        binsizes.resize(npix_unbinned,0); // This array will be clipped to npix.
        npix = 0; // Will be updated to highest non-fillvalue element in the binning table plus one.
        for (size_t ibin_spat=0 ; ibin_spat<bin_file->dim_binningtable_spat ; ibin_spat++) {
            for (size_t ibin_spec=0 ; ibin_spec<bin_file->dim_binningtable_spec ; ibin_spec++) {
                uint32_t index_cur = binning_table_read[ibin_spat*bin_file->dim_binningtable_spec+ibin_spec];
                if (index_cur != NC_FILL_UINT) {
                    if (index_cur+1 > npix) npix = index_cur+1; // Eventually update number of binned pixels.
                    binsizes[index_cur] += bin_file->prebinning_spat*bin_file->prebinning_spec; // Count bin size.
                }
                for (size_t iprebin_spat=0 ; iprebin_spat<bin_file->prebinning_spat ; iprebin_spat++) {
                    size_t ispat = ibin_spat*bin_file->prebinning_spat + iprebin_spat;
                    for (size_t iprebin_spec=0 ; iprebin_spec<bin_file->prebinning_spec ; iprebin_spec++) {
                        size_t ispec = ibin_spec*bin_file->prebinning_spec + iprebin_spec;
                        pixelpointer[ispat*ckd->dim_detector_spec+ispec] = index_cur; // Write this index, fillvalue or not.
                    }
                }
            }
        }
        binsizes.resize(npix); // Clip to correct size.
        // Construct the binned CKD with all the detector steps.
        own_binned_ckd = make_unique<CKD>(this);
        binned_ckd = own_binned_ckd.get(); // Set pointer to own binned CKD.
        // We will go to the same level as the original CKD.
        binned_ckd->lev = ckd->lev;
        // Make the CKD that we need.
        if (ckd->lev > LEVEL_DIMCAL) { // {{{
            // Copy dimensions of the detector itself.
            binned_ckd->dim_detector_spat = ckd->dim_detector_spat;
            binned_ckd->dim_detector_spec = ckd->dim_detector_spec;
            binned_ckd->dim_vp = ckd->dim_vp;
            // But save the number of pixels that belong to the binning table.
            binned_ckd->npix = npix;
            // Binned pixel mask.
            binned_ckd->mask.resize(npix,false);
            for (size_t ipix_unbinned=0 ; ipix_unbinned<npix_unbinned ; ipix_unbinned++) {
                // If a small pixel is dead, the aggregate pixel it points to is sentenced to death.
                if (ckd->mask[ipix_unbinned] && pixelpointer[ipix_unbinned] < npix) binned_ckd->mask[pixelpointer[ipix_unbinned]] = true;
            }
            binned_ckd->vp_mask = ckd->vp_mask; // Viewport mask (only relevant if we move to binned CKD after straylight).
            // Furthermore, bins with zero size are skipped.
            for (size_t ipix=0 ; ipix<npix ; ipix++) {
                if (binsizes[ipix] == 0) binned_ckd->mask[ipix] = true;
            }
            // Count number of living pixels (for convenience).
            binned_ckd->nliving = 0;
            for (size_t ipix=0 ; ipix<npix ; ipix++) if (!binned_ckd->mask[ipix]) binned_ckd->nliving++;
        } // }}}
        if (ckd->lev > LEVEL_DARKCAL) { // {{{
            // Copy skip flag.
            binned_ckd->dark_skip = ckd->dark_skip;
            // Average the relevant detector CKD if this step is not skipped.
            if (!ckd->dark_skip) {
                binned_ckd->dim_dark_order = ckd->dim_dark_order;
                handle(average_ckd(1,ckd->dim_dark_order,ckd->dark_offset,binned_ckd->dark_offset));
                handle(average_ckd(1,ckd->dim_dark_order,ckd->dark_current,binned_ckd->dark_current));
                binned_ckd->dark_nominal_temperature = ckd->dark_nominal_temperature;
            }
        } // }}}
        if (ckd->lev > LEVEL_NOISECAL) { // {{{
            // Copy skip flag.
            binned_ckd->noise_skip = ckd->noise_skip;
            // Average the relevant detector CKD if this step is not skipped.
            if (!ckd->noise_skip) {
                handle(average_ckd(1,1,ckd->noise_g,binned_ckd->noise_g));
                handle(average_ckd(1,1,ckd->noise_n,binned_ckd->noise_n));
            }
        } // }}}
        if (ckd->lev > LEVEL_NONLINCAL) { // {{{
            // Copy skip flag.
            binned_ckd->nonlin_skip = ckd->nonlin_skip;
            if (!ckd->nonlin_skip) {
                // Copy a few dimensions.
                binned_ckd->dim_nonlin_exptime = ckd->dim_nonlin_exptime;
                binned_ckd->dim_nonlin_spline = ckd->dim_nonlin_spline;
                binned_ckd->dim_nonlin_knot = ckd->dim_nonlin_knot;
                // Shape CKD-owned arrays.
                binned_ckd->nonlin_knots = ckd->nonlin_knots; // Copy the vector.
                handle(average_ckd(1,ckd->dim_nonlin_exptime*ckd->dim_nonlin_spline,ckd->nonlin_fit,binned_ckd->nonlin_fit));
                binned_ckd->nonlin_exptimes = ckd->nonlin_exptimes; // Copy the vector.
                handle(average_ckd(1,ckd->dim_nonlin_exptime,ckd->nonlin_signal_scale,binned_ckd->nonlin_signal_scale));
                binned_ckd->nonlin_order = ckd->nonlin_order;
            }
        } // }}}
        if (ckd->lev > LEVEL_PRNUCAL) { // {{{
            // Copy skip flag.
            binned_ckd->prnu_skip = ckd->prnu_skip;
            if (!ckd->prnu_skip) {
                handle(average_ckd(1,1,ckd->prnu_prnu,binned_ckd->prnu_prnu));
            }
        } // }}}
        // Straylight CKD will never be binned. Straylight is a step that happens in
        // the unbinned perspective.
        binned_ckd->stray_skip = ckd->stray_skip;

        // FOV CKD. Meanwhile, help CKD is generated with which the next steps can be simplified.
        vector<size_t> help_ispec_start;
        vector<size_t> help_ispec_end;
        size_t &nspec_total = binned_ckd->dim_fov_spec_total; // Counter of all spectral pixel together.
        nspec_total = 0;
        if (ckd->lev > LEVEL_FOVCAL) { // {{{
            // Get number of FOVs.
            binned_ckd->dim_fov = ckd->dim_fov;
            // Copy the rest of the organizational content.
            binned_ckd->fov_nfov_vp = ckd->fov_nfov_vp; // Copy the vector.
            binned_ckd->fov_act_angles = ckd->fov_act_angles; // Only needed by swathcal, which has not binning, but to prevent any chance of unexpected problems, we will copy this small vector.
            // Part of the FOV CKD is trivial for the original CKD, but non-trivial for
            // the binned CKD.
            // Shape these arrays.
            binned_ckd->fov_dims_spec.resize(ckd->dim_fov);
            binned_ckd->fov_iel_start.resize(ckd->dim_fov);
            binned_ckd->fov_ipix1.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec); // Maximum size.
            binned_ckd->fov_ipix2.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec); // Maximum size.
            binned_ckd->fov_weight1.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec); // Maximum size.
            // Help CKD.
            help_ispec_start.resize(ckd->dim_fov*ckd->dim_detector_spec); // Maximum size.
            help_ispec_end.resize(ckd->dim_fov*ckd->dim_detector_spec); // Maximum size.

            // Make pointers to target CKD.
            size_t *ipix1_cur = binned_ckd->fov_ipix1.data();
            size_t *ipix2_cur = binned_ckd->fov_ipix2.data();
            double *weight1_cur = binned_ckd->fov_weight1.data();
            size_t *help_ispec_start_cur = help_ispec_start.data();
            size_t *help_ispec_end_cur = help_ispec_end.data();
            // Pointer to source CKD.
            double *fov_ispat_cur = ckd->fov_ispat.data();
            for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
                size_t sz_pol_m; // To verify that the spectra for S+ and S- are equally large.
                // In fact, they should also have the same spectral binning, although that is
                // algorithmically not necessary. The same spectral binning would mean that S+ and
                // S- both approximately the same wavelength sampling, but it is not really the same,
                // because of the wavelength CKD. Having larger differences in the wavelength sampling
                // may be harmful due to interpolation errors on target wavelengths when using the
                // antisymmetric demodulation. If non-antisymmetric demodulation is chosen, the spectra
                // need not to be equally large. Currently, the Spectra structure likes having same-sized
                // spectra for S+ and S- for the same FOV.
                
                for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                    // The idea is to calculate the ispats of the original CKD.
                    // Hidden for-loop. I want to execute a part in the beginning again after
                    // the last iteration.
                    size_t ispec_orig = 0;
                    size_t ispec_bin_start = NC_FILL_UINT64;
                    size_t ispec_bin_end = NC_FILL_UINT64;
                    size_t binsize = 0;
                    double agg_ispat = 0.0;
                    size_t sz = 0;

                    while (true) { // Full iterators from 0 to dim_detector_spec and a partial one at dim_detector_spec.

                        // Step one. Close off previous cell. {{{
                        // This is when your index grows to equal to the ending cell.
                        // Then, the last cell is closed when ispec_orig
                        // reaches ckd->dim_detector_spec.
                        if (ispec_orig == ispec_bin_end) {

                            // If something illegal happens, an error is raised. Filling fillvalues would be an option,
                            // but it can become difficult to follow then.
                            double ispat_avg = agg_ispat / binsize;
                            check_error(ispat_avg < 0.0 || ispat_avg > static_cast<double>(ckd->dim_detector_spat)-1.0,"Program error: Spectrum outside detector. Possibly binning table and CKD do not match. This should have been caught before.");
                            size_t ispat_rounded = static_cast<size_t>(ispat_avg+0.5);
                            size_t ipix_main; // Output of explore: Main binned pixel index.
                            size_t ispat_start; // Output of explore, starting spatial index.
                            size_t ispat_end; // Output of explore, ending spatial index.
                            handle(explore_cell(ispat_rounded,ispec_bin_start,ispec_bin_end,ipix_main,ispat_start,ispat_end));
                            double ispat_center_main = 0.5*static_cast<double>(ispat_start+ispat_end-1);
                            if (ispat_avg < ispat_center_main) {
                                // Desire to move down.
                                if (ispat_start == 0) {
                                    // Impossible to move down. Take just the main cell.
                                    *ipix1_cur = ipix_main;
                                    *ipix2_cur = ipix_main;
                                    *weight1_cur = 1.0;
                                } else {
                                    // Move down.
                                    size_t ipix_secondary;
                                    size_t ispat_start_secondary;
                                    size_t ispat_end_secondary;
                                    handle(explore_cell(ispat_start-1,ispec_bin_start,ispec_bin_end,ipix_secondary,ispat_start_secondary,ispat_end_secondary));
                                    check_error(ispat_end_secondary != ispat_start,"Program error: It looks like there a wall with only one side.");
                                    double ispat_center_secondary = 0.5*static_cast<double>(ispat_start_secondary+ispat_end_secondary-1);
                                    *ipix1_cur = ipix_secondary; // Lower.
                                    *ipix2_cur = ipix_main; // Higher.
                                    *weight1_cur = (ispat_center_main-ispat_avg) / (ispat_center_main-ispat_center_secondary);
                                }
                            } else {
                                // Desire to move up.
                                if (ispat_end == ckd->dim_detector_spat) {
                                    // Impossible to move up. Take just the main cell.
                                    *ipix1_cur = ipix_main;
                                    *ipix2_cur = ipix_main;
                                    *weight1_cur = 1.0;
                                } else {
                                    // Move up.
                                    size_t ipix_secondary;
                                    size_t ispat_start_secondary;
                                    size_t ispat_end_secondary;
                                    handle(explore_cell(ispat_end,ispec_bin_start,ispec_bin_end,ipix_secondary,ispat_start_secondary,ispat_end_secondary));
                                    check_error(ispat_start_secondary != ispat_end,"Program error: It looks like there a wall with only one side.");
                                    double ispat_center_secondary = 0.5*static_cast<double>(ispat_start_secondary+ispat_end_secondary-1);
                                    *ipix1_cur = ipix_main; // Lower.
                                    *ipix2_cur = ipix_secondary; // Higher.
                                    *weight1_cur = (ispat_center_secondary-ispat_avg) / (ispat_center_secondary-ispat_center_main);
                                }
                            }

                            if (ipol == 0) {
                                // Write cell borders to help CKD.
                                *help_ispec_start_cur = ispec_bin_start;
                                *help_ispec_end_cur = ispec_bin_end;
                            } else {
                                // Verify consistency between the two polarization directions.
                                check_error(ispec_bin_start != *help_ispec_start_cur,"Error: Spectral binning seems to differ between S+ and S- for FOV %zu. S- starts a cell at index %zu. S+ starts the same cell at index %zu",ifov,*help_ispec_start_cur,ispec_bin_start);
                                check_error(ispec_bin_end != *help_ispec_end_cur,"Error: Spectral binning seems to differ between S+ and S- for FOV %zu. S- end a cell at index %zu. S+ end the same cell at index %zu",ifov,*help_ispec_end_cur,ispec_bin_end);
                            }

                            sz++;
                            ipix1_cur++;
                            ipix2_cur++;
                            weight1_cur++;
                            help_ispec_start_cur++;
                            help_ispec_end_cur++;
                        } // }}}
                        // Step two. Stop if your iterator is no longer on the detector.
                        if (ispec_orig == ckd->dim_detector_spec) break;
                        // Step three. Open a new cell. // {{{
                        // Nice to know always.
                        if (ispec_orig == 0 || ispec_orig == ispec_bin_end) {
                            if (ispec_orig == 0) {
                                ispec_bin_start = 0;
                                ispec_bin_end = 0; // This is an exploration starting point.
                            } else {
                                ispec_bin_start = ispec_bin_end;
                            }
                            // Reset counter.
                            binsize = 0;
                            agg_ispat = 0.0;
                        } // }}}

                        // Step four. Explore current cell.
                        // If the cell is just opened, the size will be defined.
                        // If the cell is already open, the explore routine will assert
                        // that the cell stays with the same size.
                        size_t ipix_ignore;
                        size_t ispat_ignore1;
                        size_t ispat_ignore2;
                        size_t ispat = static_cast<size_t>(fov_ispat_cur[ispec_orig]+0.5);
                        handle(explore_cell(ispat,ispec_bin_start,ispec_bin_end,ipix_ignore,ispat_ignore1,ispat_ignore2));
                        binsize++;
                        // Step five. Count current cell.
                        agg_ispat += fov_ispat_cur[ispec_orig];
                        ispec_orig++; // Progress iterator.
                    } // Clumsy while-true for-loop.
                    // Still inside the loop over polarizations.
                    if (ipol == 0) sz_pol_m = sz;
                    else {
                        check_error(sz != sz_pol_m,"Error: S+ and S- spectra are not equally long for FOV %zu. S- has size %zu and S+ has size %zu",ifov,sz_pol_m,sz);
                    }

                    // Rewind pointers that do not have the polarization dimension after
                    // the first polarization loop iteration..
                    if (ipol != DIM_POL-1) {
                        help_ispec_start_cur -= sz;
                        help_ispec_end_cur -= sz;
                    }

                    // Increase source CKD pointer.
                    fov_ispat_cur += ckd->dim_detector_spec;
                }
                // Outside the pol loop, still in the fov loop.
                binned_ckd->fov_iel_start[ifov] = DIM_POL*nspec_total; // Before adding current size.
                binned_ckd->fov_dims_spec[ifov] = sz_pol_m; // That is the remaining counter.
                nspec_total += sz_pol_m; // Count up currentsize.
            } // Loop over FOVs.

            // Encapsulate target arrays.
            // nspec_total is now the fold of the spectral dimension and the FOV dimension, but not the POL dimension.
            binned_ckd->fov_ipix1.resize(DIM_POL*nspec_total);
            binned_ckd->fov_ipix2.resize(DIM_POL*nspec_total);
            binned_ckd->fov_weight1.resize(DIM_POL*nspec_total);
            help_ispec_start.resize(nspec_total);
            help_ispec_end.resize(nspec_total);

        } // }}}

        // Swath definition.
        if (ckd->lev > LEVEL_SWATHCAL) {
            // Copy the information.
            binned_ckd->swath_skip = ckd->swath_skip;
            if (!ckd->swath_skip) {
                binned_ckd->swath_swathvectors = ckd->swath_swathvectors; // Copy the entire contents.
                binned_ckd->swath_vectorplane_normals = ckd->swath_vectorplane_normals; // Copy the entire contents.
            }
        }

        // Wavelength calibration.
        if (ckd->lev > LEVEL_WAVECAL) { // {{{
            binned_ckd->wave_spectra.resize(nspec_total*DIM_POL);
            binned_ckd->wave_target.resize(nspec_total);
            // Average the wavelengths.
            size_t *ispec_start = help_ispec_start.data();
            size_t *ispec_end = help_ispec_end.data();
            double *wave_spectra_out = binned_ckd->wave_spectra.data();
            double *wave_target_out = binned_ckd->wave_target.data();
            double *wave_spectra_in = ckd->wave_spectra.data();
            for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
                size_t &nbin = binned_ckd->fov_dims_spec[ifov];
                for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                    for (size_t ibin=0 ; ibin<nbin ; ibin++) {
                        double agg = 0.0;
                        for (size_t ispec=ispec_start[ibin] ; ispec<ispec_end[ibin] ; ispec++) {
                            agg += wave_spectra_in[ipol*ckd->dim_detector_spec+ispec];
                        }
                        wave_spectra_out[ipol*nbin+ibin] = agg / static_cast<double>(ispec_end[ibin]-ispec_start[ibin]);
                    }
                }
                for (size_t ibin=0 ; ibin<nbin ; ibin++) {
                    wave_target_out[ibin] = 0.5*(wave_spectra_out[ibin] + wave_spectra_out[nbin+ibin]);
                }
                // Progress all pointers.
                ispec_start += nbin;
                ispec_end += nbin;
                wave_spectra_out += DIM_POL*nbin;
                wave_target_out += nbin;
                wave_spectra_in += DIM_POL*ckd->dim_detector_spec;
            }
        } // }}}

        // Radiometric calibration.
        if (ckd->lev > LEVEL_RADCAL) { // {{{
            binned_ckd->rad_spectra.resize(nspec_total*DIM_POL);
            binned_ckd->rad_skip = ckd->rad_skip;
            // I think that the harmonic mean is best.
            size_t *ispec_start = help_ispec_start.data();
            size_t *ispec_end = help_ispec_end.data();
            double *rad_spectra_out = binned_ckd->rad_spectra.data();
            double *rad_spectra_in = ckd->rad_spectra.data();
            for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
                size_t &nbin = binned_ckd->fov_dims_spec[ifov];
                for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                    for (size_t ibin=0 ; ibin<nbin ; ibin++) {
                        double agg = 0.0;
                        for (size_t ispec=ispec_start[ibin] ; ispec<ispec_end[ibin] ; ispec++) {
                            agg += 1.0/rad_spectra_in[ipol*ckd->dim_detector_spec+ispec];
                        }
                        rad_spectra_out[ipol*nbin+ibin] = static_cast<double>(ispec_end[ibin]-ispec_start[ibin]) / agg;
                        // Progress all pointers.
                    }
                }
                ispec_start += nbin;
                ispec_end += nbin;
                rad_spectra_out += DIM_POL*nbin;
                rad_spectra_in += DIM_POL*ckd->dim_detector_spec;
            }
        } // }}}
    }
    return 0;
}

int Binningtable::average_ckd( // {{{
    size_t dim_slower, // Product of dimensions slower than pixels.
    size_t dim_quicker, // Product of dimensions quicker than pixels.
    vector<double> &src, // Source vector from original CKD.
    vector<double> &dest // Destination vector from binned CKD.
)
{
    dest.resize(dim_slower*npix*dim_quicker,0.0);
    check_error(src.size() != dim_slower*npix_unbinned*dim_quicker,"Error: Source array has not the correct size: Expected %zu, Actual %zu.",dim_slower*npix_unbinned*dim_quicker,src.size());

    for (size_t ipix_unbinned=0 ; ipix_unbinned<npix_unbinned ; ipix_unbinned++) {
        uint32_t &ipix = pixelpointer[ipix_unbinned];
        if (ipix < npix) {
            for (size_t islow=0 ; islow<dim_slower ; islow++) {
                for (size_t iquick=0 ; iquick<dim_quicker ; iquick++) {
                    dest[islow*npix*dim_quicker+ipix*dim_quicker+iquick] += src[islow*npix_unbinned*dim_quicker+ipix_unbinned*dim_quicker+iquick];
                }
            }
        }
    }
    for (size_t ipix=0 ; ipix<npix ; ipix++) {
        if (binsizes[ipix] == 0) {
            // Prevent zero-division.
            for (size_t islow=0 ; islow<dim_slower ; islow++) {
                for (size_t iquick=0 ; iquick<dim_quicker ; iquick++) {
                    dest[islow*npix*dim_quicker+ipix*dim_quicker+iquick] = NC_FILL_DOUBLE;
                }
            }
        } else {
            for (size_t islow=0 ; islow<dim_slower ; islow++) {
                for (size_t iquick=0 ; iquick<dim_quicker ; iquick++) {
                    dest[islow*npix*dim_quicker+ipix*dim_quicker+iquick] /= binsizes[ipix];
                }
            }
        }
    }

    return 0;

} // }}}

int Binningtable::explore_cell( // {{{
    size_t ispat,
    size_t ispec_start,
    size_t &ispec_end,
    size_t &ipix,
    size_t &ispat_start,
    size_t &ispat_end
)
{

    // Walk at ispat from left to right inside the cell and feel the
    // walls on the left and the right.
    ipix = pixelpointer[ispat*binned_ckd->dim_detector_spec+ispec_start]; // On the left of the cell. That must be your cell.
    check_error(ipix >= npix,"Error: Spectra goes through thrown-away pixels.");
    // Feel the wall on the left.
    if (ispec_start != 0) { // If ispec_start is zero, there is an outside wall left of you.
        check_error(pixelpointer[ispat*binned_ckd->dim_detector_spec+ispec_start-1] == ipix,"Error: Left wall of cell not on the expected place."); // No wall where one is expected.
    }
    // If the ending pixel was equal to the starting pixel, explore the size in spectral direction instead.
    if (ispec_end == ispec_start) {
        // Explore.
        for (ispec_end=ispec_start+1 ; ispec_end<=binned_ckd->dim_detector_spec ; ispec_end++) {
            if (ispec_end == binned_ckd->dim_detector_spec) break; // Outside wall.
            if (pixelpointer[ispat*binned_ckd->dim_detector_spec+ispec_end] != ipix) break; // Hit an interior wall.
        }
    } else {
        // A spectral domain is given, so now, check if this is correct.
        // Walk through the cell hoping to keep your nose (not hitting an unexpected wall).
        for (size_t ispec_explore=ispec_start+1 ; ispec_explore<ispec_end ; ispec_explore++) {
            check_error(pixelpointer[ispat*binned_ckd->dim_detector_spec+ispec_explore] != ipix,"Error: Unexpected right wall in cell."); // Hit the wall.
        }
        // Feel the wall on the right.
        if (ispec_end != binned_ckd->dim_detector_spec) { // If ispec_bin_end is dim_detector_spec, there is an outside wall right of you.
            check_error(pixelpointer[ispat*binned_ckd->dim_detector_spec+ispec_end] == ipix,"Error: Right wall was expected and is not there."); // No wall where one is expected.
        }
    }
    // Explore vertically to get the average spatial index.
    // We assume that the cell is rectangular. Otherwise, we have a bigger problem.
    ispat_start = 0; // When not broken, the cell starts at the bottom.
    ispat_end = binned_ckd->dim_detector_spat; // When not broken, the cell ends at the top.
    for (size_t ispat_cur=ispat ; ispat_cur > 0 ; ispat_cur--) {
        // We avoided using a -1 with unsigned integers.
        size_t ispat_explore = ispat_cur-1;
        if (pixelpointer[ispat_explore*binned_ckd->dim_detector_spec+ispec_start] != ipix) {
            // Hit the wall, recover your nose and that is the start.
            ispat_start = ispat_cur;
            break;
        }
    }
    for (size_t ispat_explore=ispat+1 ; ispat_explore<binned_ckd->dim_detector_spat ; ispat_explore++) {
        if (pixelpointer[ispat_explore*binned_ckd->dim_detector_spec+ispec_start] != ipix) {
            // Hit the wall, keep your bleeding nose and that is the end.
            // This is because we always define the end as the first excluded pixel, so therefore, your nose should stay bleeding here.
            ispat_end = ispat_explore;
            break;
        }
    }

    return 0;

} // }}}

