// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "array.h"
#include "matrix.h"
#include "lininv.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "nonlincal.h"

// Settings functions.
Settings_nonlincal::Settings_nonlincal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "nonlin";
    output_level_max = 1;
} // }}}
Settings_nonlincal::~Settings_nonlincal() {}
int Settings_nonlincal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(nfrac); // Number of fraction to break up L1A images.
    recognize_setting(lmin); // Lower border (in measured signal) of linear domain.
    recognize_setting(lmax); // Upper border (in measured signal) of linear domain.
    recognize_setting_vector(exposure_time); // Exposure times for which the non-linearity is evaluated.
    recognize_setting(exptime_tol); // Relative tolerance for exposure time grouping.
    recognize_setting(scale_median); // Flag to scale the spline knots with the median of the saturated image. Set false to use individual-pixel saturation levels.
    recognize_setting_vector(knots); // B-spline knots in signal domain.
    recognize_setting(order); // B-spline order.
    recognize_setting(mask_lin_chi2_max); // Maximum allowed residual chi squared on the linear fit during non-linearity.
    recognize_setting(mask_chi2_max); // Maximum allowed residual chi squared on function fit of the non-linearity correction.

    return 0;

} // }}}

// Constructor for the Nonlincal CKD structure.
Nonlincal::Nonlincal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("nonlin");
    set = make_unique<Settings_nonlincal>(this);
    Processor::set = set.get();
    own_skip = &ckd->nonlin_skip;
} // }}}
Nonlincal::~Nonlincal() {}

int Nonlincal::process_init( // {{{
)
{

    // Verify if the settings are legal.
    check_error(set->lmin == NC_FILL_DOUBLE,"Error: Missing minimum signal in linear region. Setting 'lmin'");
    check_error(set->lmax == NC_FILL_DOUBLE,"Error: Missing maximum signal in linear region. Setting 'lmax'");
    check_error(set->exposure_time.size() == 0,"Error: No operational exposure time given. Setting 'exposure_time'");
    check_error(set->exptime_tol < 1.0,"Error: Exposure time ratio tolerance should be at least one. Setting 'exptime_tol'");
    check_error(set->knots.size() < 2,"Error: B-splines need at least two knots. Setting 'knots'");
    for (size_t iknot=1 ; iknot<set->knots.size() ; iknot++) {
        check_error(set->knots[iknot] < set->knots[iknot-1],"Error: Knots must be set in ascending order.");
    }
    check_error(set->order < 1,"Error: B-spline degrees of freedom must be at least one, for a constant correction per knot interval. Setting 'order'");

    // Verify existence of necessary GSE data.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        check_error(l1a_instances_all[il1a]->illumination_level == NC_FILL_DOUBLE,"Error: Missing illumination level as GSE-data from '%s'.",l1a_instances_all[il1a]->filename.c_str());
    }

    // Create the fitting B-spline, meanwhile leaving all essenntial information
    // on the 'nonlin'-prefixed variables in the CKD structure.
    ckd->dim_nonlin_knot = set->knots.size();
    ckd->nonlin_knots.resize(ckd->dim_nonlin_knot);
    memcpy(ckd->nonlin_knots.data(),set->knots.data(),ckd->dim_nonlin_knot*sizeof(double));
    ckd->nonlin_order = set->order;
    b = make_unique<Bspline>(this,ckd->nonlin_order,ckd->dim_nonlin_knot,ckd->nonlin_knots.data());
    ckd->dim_nonlin_spline = b->nspline;

    // Recycled element: The abscissa.
    abscissa.resize(nl1a_total);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) abscissa[il1a] = l1a_instances_all[il1a]->exposure_time*l1a_instances_all[il1a]->illumination_level;

    // Find out how many exposure times are used.
    ckd->dim_nonlin_exptime = set->exposure_time.size();
    ckd->nonlin_exptimes.resize(ckd->dim_nonlin_exptime);
    memcpy(ckd->nonlin_exptimes.data(),set->exposure_time.data(),ckd->dim_nonlin_exptime*sizeof(double));

    // Verify that no exposure time is inside the tolerance limit of two different selected exposure times.
    double tol2 = pow(set->exptime_tol,2.0);
    for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
        double &exptime1 = ckd->nonlin_exptimes[iexptime];
        for (size_t iexptime_check=0 ; iexptime_check<iexptime ; iexptime_check++) {
            double &exptime0 = ckd->nonlin_exptimes[iexptime_check];
            check_error(exptime1 / exptime0 <= tol2 && exptime0 / exptime1 <= tol2,"Settings error: Tolerance areas of two selected exposure times overlap: %.7f and %.7f.",exptime0,exptime1);
        }
    }

    // First, sort out which files belong to a chosen exposure time (iexptime). If it not not
    // associated to an exposure time, mark dim_nonlin_exptime.
    // Meanwhile, count the number of files per exposure time.
    vector<size_t> iexptime_l1a(nl1a_total);
    vector<size_t> nl1a_exptime(ckd->dim_nonlin_exptime,0);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        double &exptime_l1a = l1a_instances_all[il1a]->exposure_time;
        bool found = false;
        for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
            double &exptime_select = ckd->nonlin_exptimes[iexptime];
            if (exptime_l1a / exptime_select <= set->exptime_tol && exptime_select / exptime_l1a <= set->exptime_tol) {
                found = true;
                iexptime_l1a[il1a] = iexptime;
                nl1a_exptime[iexptime]++;
                break;
            }
        }
        if (!found) iexptime_l1a[il1a] = ckd->dim_nonlin_exptime;
    }

    // Verify that every exposure time has emough images.
    for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
        check_error(nl1a_exptime[iexptime] < ckd->dim_nonlin_spline,"Error: Unable to perform non-linearity fit for expousre tim %.3f, because there are only %zu measurements and %zu degrees of freedom in the fit.",ckd->nonlin_exptimes[iexptime],nl1a_exptime[iexptime],ckd->dim_nonlin_spline);
    }

    // The following batches will be formed.
    // - Highest illumination levels of chosen exposure times.
    // --- dim_nonlin_exptime batches, not chunked.
    // - All files.
    // --- One batch, chunked.
    // - All files for chosen exposure time (dim_nonlin_exptime batches of multiple files).
    // --- dim_nonlin_exptime batches, chunked.

    // Add administration to the batches: The phase that belongs to the batch and the exposure time
    // that belongs to the batch (or dim_nonlin_exptime if it is not specific to an exposure time).

    // Total number of batches: dim_nonlin_exptime + nfrac + dim_nonlin_exptime*nfrac
    nbatch = ckd->dim_nonlin_exptime + set->nfrac + ckd->dim_nonlin_exptime*set->nfrac;
    batches.resize(nbatch);
    phase_batches.resize(nbatch);
    iexptime_batches.resize(nbatch);
    ifrac_batches.resize(nbatch);
    // Reverse this administration for here.
    vector<size_t> ibatch_medianphase(ckd->dim_nonlin_exptime);
    vector<size_t> ibatch_linphase(set->nfrac);
    vector<size_t> ibatch_nonlinphase(set->nfrac*ckd->dim_nonlin_exptime); // Slow dimension is fraction. Quick dimension is exposure time.

    // First dim_nonlin_exptime batches are the whole images for scaling.
    size_t ibatch = 0; // Counter.
    // Of course, here ibatch is iexptime, but it is better to keep it clear which index is meant.
    for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
        ibatch_medianphase[iexptime] = ibatch;
        phase_batches[ibatch] = phase_median;
        iexptime_batches[ibatch] = iexptime;
        ifrac_batches[ibatch] = set->nfrac; // Not used.
        ibatch++;
    }
    // Next there are nfrac series of lin, nonlin, nonlin ..., where number of nonlins is dim_nonlin_exptime.
    for (size_t ifrac=0 ; ifrac<set->nfrac ; ifrac++) {
        ibatch_linphase[ifrac] = ibatch;
        phase_batches[ibatch] = phase_lin;
        iexptime_batches[ibatch] = ckd->dim_nonlin_exptime; // Not usd.
        ifrac_batches[ibatch] = ifrac;
        ibatch++;
        for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
            ibatch_nonlinphase[ifrac*ckd->dim_nonlin_exptime+iexptime] = ibatch;
            phase_batches[ibatch] = phase_nonlin;
            iexptime_batches[ibatch] = iexptime;
            ifrac_batches[ibatch] = ifrac;
            ibatch++;
        }
    }

    // Set the sizes of the batches.
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        Batch &bat = batches[ibatch];
        if (phase_batches[ibatch] == phase_median) {
            bat.setSize(1);
            // Only the L1A* has to be set. The rest stays at default.
            // This is done later.
            bat.remove[0] = false;
        } else if (phase_batches[ibatch] == phase_lin) {
            bat.setSize(nl1a_total);
            bat.ipix_start = ckd->npix * ifrac_batches[ibatch] / set->nfrac;
            bat.ipix_end = ckd->npix * (ifrac_batches[ibatch]+1) / set->nfrac;
            for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
                bat.l1a[il1a] = l1a_instances_all[il1a].get();
                bat.remove[il1a] = iexptime_l1a[il1a] == ckd->dim_nonlin_exptime;
                // Add remove flag (true for everything except for those that have an exposure time).
            }
        } else if (phase_batches[ibatch] == phase_nonlin) {
            size_t &iexptime = iexptime_batches[ibatch];
            size_t &nl1a_cur = nl1a_exptime[iexptime];
            bat.setSize(nl1a_cur);
            bat.ipix_start = ckd->npix * ifrac_batches[ibatch] / set->nfrac;
            bat.ipix_end = ckd->npix * (ifrac_batches[ibatch]+1) / set->nfrac;
            for (size_t il1a=0 ; il1a<nl1a_cur ; il1a++) {
                // L1A* will be set when browsing through L1A instances.
                bat.remove[il1a] = true;
            }
        }
    }

    // Now, store the L1A pointers in the nonlin and the median phases.
    vector<double> illum_max(ckd->dim_nonlin_exptime,0.0);
    vector<size_t> cnt_exptime(ckd->dim_nonlin_exptime,0); // Batch element counter.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        L1A *l1a_cur = l1a_instances_all[il1a].get();
        size_t &iexptime = iexptime_l1a[il1a];
        if (iexptime != ckd->dim_nonlin_exptime) {
            // Append to nonlin list.
            size_t &cnt = cnt_exptime[iexptime];
            for (size_t ifrac=0 ; ifrac<set->nfrac ; ifrac++) {
                size_t &ibatch = ibatch_nonlinphase[ifrac*ckd->dim_nonlin_exptime+iexptime];
                Batch &bat = batches[ibatch];
                bat.l1a[cnt] = l1a_cur;
            }
            cnt++;
            // If this is the largest illumination level in the exposure time, overwrite the L1A* in
            // the median phase.
            if (l1a_cur->illumination_level > illum_max[iexptime]) {
                illum_max[iexptime] = l1a_cur->illumination_level;
                // Overwrite the L1A* (or set it if this is the first time.
                size_t &ibatch = ibatch_medianphase[iexptime];
                Batch &bat = batches[ibatch];
                bat.l1a[0] = l1a_cur;
            }
        }
    }

    // Shape the rest of the CKD.
    ckd->nonlin_signal_scale.resize(ckd->npix*ckd->dim_nonlin_exptime);
    nonlin_lin_slope.resize(ckd->npix);
    nonlin_lin_chi2.resize(ckd->npix);
    ckd->nonlin_fit.resize(ckd->npix*ckd->dim_nonlin_exptime*ckd->dim_nonlin_spline);
    nonlin_chi2.resize(ckd->npix*ckd->dim_nonlin_exptime);
    // Detailed output.
    if (set->output_level >= 1) {
        detailed_output_series.resize(ckd->dim_nonlin_exptime);
        for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
            DetailedOutputNonlin &det_out = detailed_output_series[iexptime];
            size_t &nl1a_cur = nl1a_exptime[iexptime];
            det_out.nl1a = nl1a_cur;
            det_out.signal_corr.resize(ckd->npix*nl1a_cur);
            det_out.signal_meas.resize(ckd->npix*nl1a_cur);
        }
    }

    // Set pointers to diagnostic CKD (owned by processor).
    ckd->diag_nonlin_lin_slope = nonlin_lin_slope.data();
    ckd->diag_nonlin_lin_chi2 = nonlin_lin_chi2.data();
    ckd->diag_nonlin_chi2 = nonlin_chi2.data();

    return 0;

} // }}}

// Non-linearity correction calibration software.
// This is the planned protocol.
int Nonlincal::process_batch(size_t ibatch, const Calibration_options& opt)
{

    Batch &bat = batches[ibatch];
    if (phase_batches[ibatch] == phase_median) {
        handle(process_median(iexptime_batches[ibatch])); // No fraction.
    } else if (phase_batches[ibatch] == phase_lin) {
        handle(process_lin(bat.ipix_start,bat.ipix_end)); // No exptime.
    } else if (phase_batches[ibatch] == phase_nonlin) {
        handle(process_nonlin(bat.ipix_start,bat.ipix_end,iexptime_batches[ibatch]));
    } else {
        raise_error("Program error: Invalid non-linearity phase.");
    }

    return 0;

} // }}}

int Nonlincal::process_median( // {{{
    size_t iexptime
)
{

    L1A *l1a = l1a_instances[0]; // There is only one file.
    if (set->scale_median) {
        double med;
        median(ckd->npix,l1a->image,med);
        for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) ckd->nonlin_signal_scale[ipix*ckd->dim_nonlin_exptime+iexptime] = med;
    } else {
        for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) ckd->nonlin_signal_scale[ipix*ckd->dim_nonlin_exptime+iexptime] = l1a->image[ipix];
    }

    return 0;

} // }}}

int Nonlincal::process_lin( // {{{
    size_t ipix_start,
    size_t ipix_end
)
{

    // Linear result pointers.
    double *linslopeptr = &nonlin_lin_slope[ipix_start];
    double *linchi2ptr = &nonlin_lin_chi2[ipix_start];

    percentagelog_open("Calculating linear fit per pixel");
    for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {

        percentagelog_progress(ipix,ckd->npix);

        // Put fillvalues everywhere.
        *linslopeptr = NC_FILL_DOUBLE;
        *linchi2ptr = NC_FILL_DOUBLE;

        // Breakable do-once loop.
        for (bool ctd=true ; ctd ; ctd=false) {

            // Die if you are already dead.
            if (ckd->mask[ipix]) break;
            ckd->mask[ipix] = true; // Unless you survive until the end of the loop.

            // Gather the measurements at this pixel index.
            vector<double> signal(nl1a);
            vector<double> noise(nl1a);
            vector<bool> mask(nl1a);
            for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                signal[il1a] = l1a_instances[il1a]->image[ipix];
                noise[il1a] = l1a_instances[il1a]->noise[ipix];
                mask[il1a] = l1a_instances[il1a]->pixelmask[ipix];
            }
            // Filter for signals between the chosen range.
            size_t nl1a_lin = 0;
            size_t nl1a_lin_healthy = 0;
            vector<double> abscissa_lin(nl1a); // Maximum allocation size.
            vector<double> signal_lin(nl1a); // Maximum allocation size.
            vector<double> s_y_lin(nl1a); // Maximum allocation size.
            vector<bool> mask_lin(nl1a); // Maximum allocation size.
            for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                double &sig = signal[il1a];
                if (sig > set->lmin && sig < set->lmax) {
                    abscissa_lin[nl1a_lin] = abscissa[il1a];
                    signal_lin[nl1a_lin] = sig;
                    s_y_lin[nl1a_lin] = pow(noise[il1a],2.0);
                    mask_lin[nl1a_lin] = mask[il1a];
                    if (!mask_lin[nl1a_lin]) nl1a_lin_healthy++;
                    nl1a_lin++;
                }
            }

            // Encapsulate vectors to their intended sizes.
            abscissa_lin.resize(nl1a_lin);
            signal_lin.resize(nl1a_lin);
            s_y_lin.resize(nl1a_lin);
            mask_lin.resize(nl1a_lin);

            // Nothing is possible if there are no reliable measurments.
            if (nl1a_lin_healthy < 1) break;

            // Perform the linear fit. Note that if it fails, no mess is created in the result before
            // returning.
            if (linear_invert(1,nl1a_lin,abscissa_lin.data(),OPT_DIAG,s_y_lin.data(),&mask_lin,signal_lin.data(),linslopeptr) != 0) break;
            *linchi2ptr = 0.0;
            for (size_t il1a=0 ; il1a<nl1a_lin ; il1a++) {
                if (mask_lin[il1a]) continue;
                *linchi2ptr += pow((abscissa_lin[il1a] * *linslopeptr - signal_lin[il1a]),2.0)/s_y_lin[il1a] / (nl1a_lin_healthy-1);
            }

            // Apply mask criterion: Linear chi squared.
            if (set->mask_lin_chi2_max != NC_FILL_DOUBLE && *linchi2ptr > set->mask_lin_chi2_max) break;

            // Survived to the end of the loop.
            ckd->mask[ipix] = false; // Unless you survive until the end of the loop.

        }

        // Progress pointers.
        linslopeptr++;
        linchi2ptr++;

    }
    int percentage_end = 100*ipix_end/ckd->npix;
    percentagelog_close(percentage_end);

    return 0;

} // }}}

int Nonlincal::process_nonlin( // {{{
    size_t ipix_start,
    size_t ipix_end,
    size_t iexptime
)
{

    // Pointers to earlier calculated results.
    double *linslopeptr = &nonlin_lin_slope[ipix_start];
    double *scaleptr = &ckd->nonlin_signal_scale[ipix_start*ckd->dim_nonlin_exptime+iexptime];
    // Pointers to results.
    double *resptr = &ckd->nonlin_fit[ipix_start*ckd->dim_nonlin_exptime*ckd->dim_nonlin_spline+iexptime*ckd->dim_nonlin_spline];
    double *chi2ptr = &nonlin_chi2[ipix_start*ckd->dim_nonlin_exptime+iexptime];

    percentagelog_open("Calculating non-linearity fit per pixel");
    for (size_t ipix=ipix_start ; ipix<ipix_end ; ipix++) {

        percentagelog_progress(iexptime*ipix_end+ipix+(ckd->dim_nonlin_exptime-iexptime-1)*ipix_start,ckd->npix*ckd->dim_nonlin_exptime);

        // Fill fillvalues everywhere.
        for (size_t ispline=0 ; ispline<ckd->dim_nonlin_spline ; ispline++) resptr[ispline] = NC_FILL_DOUBLE;
        *chi2ptr = NC_FILL_DOUBLE;
        if (set->output_level >= 1) {
            DetailedOutputNonlin &det_out = detailed_output_series[iexptime];
            for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                det_out.signal_meas[ipix*nl1a+il1a] = NC_FILL_DOUBLE;
                det_out.signal_corr[ipix*nl1a+il1a] = NC_FILL_DOUBLE;
            }
        }

        // Breakable do-once loop.
        for (bool ctd=true ; ctd ; ctd=false) {

            // Die if you are already dead.
            if (ckd->mask[ipix]) break;
            ckd->mask[ipix] = true; // Unless you survive until the end of the loop.

            vector<double> signal_corrected_scaled(nl1a);
            vector<double> signal_measured(nl1a);
            vector<double> s_y_measured(nl1a);
            vector<bool> mask_measured(nl1a);
            size_t nl1a_healthy = 0;
            for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                L1A *l1a = l1a_instances[il1a];
                signal_corrected_scaled[il1a] = *linslopeptr * l1a->exposure_time * l1a->illumination_level / (*scaleptr);
                signal_measured[il1a] = l1a->image[ipix];
                s_y_measured[il1a] = pow(l1a->noise[ipix],2.0);
                mask_measured[il1a] = l1a->pixelmask[ipix];
                if (!mask_measured[il1a]) nl1a_healthy++;
            }
            if (set->output_level >= 1) {
                // Write measured and unscaled corrected signal.
                DetailedOutputNonlin &det_out = detailed_output_series[iexptime];
                for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                    det_out.signal_meas[ipix*nl1a+il1a] = signal_measured[il1a];
                    det_out.signal_corr[ipix*nl1a+il1a] = signal_corrected_scaled[il1a]*(*scaleptr);
                }
            }

            // Perform the fit.
            vector<double> jac(ckd->dim_nonlin_spline*nl1a);
            b->jaccalc(nl1a,signal_corrected_scaled.data(),jac.data());
            // Add hard constraint at zero. Zero signal should mean zero mismatch.
            double zero = 0.0; // Used as abscissa value and as measurement.
            vector<double> jac_zero(ckd->dim_nonlin_spline);
            b->jaccalc(1,&zero,jac_zero.data());
            if (constrained_linear_invert(ckd->dim_nonlin_spline,nl1a,1,jac.data(),jac_zero.data(),OPT_DIAG,s_y_measured.data(),&mask_measured,signal_measured.data(),&zero,resptr) != 0) break;
            *chi2ptr = 0.0;
            vector<double> mod(nl1a);
            Matrix::matmul_rl_fold_slow(ckd->dim_nonlin_spline,nl1a,jac.data(),resptr,mod.data());
            for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
                if (mask_measured[il1a]) continue;
                *chi2ptr += pow(mod[il1a]-signal_measured[il1a],2.0)/s_y_measured[il1a] / (nl1a_healthy-ckd->dim_nonlin_spline+1); // The plus one is the hard constraint, reducing the DFS by one (or adding a measurement with no bias).
            }

            // Apply mask criterion: Nonlinear chi squared.
            if (set->mask_chi2_max != NC_FILL_DOUBLE && *chi2ptr > set->mask_chi2_max) break;

            // Survived to the end of the loop.
            ckd->mask[ipix] = false;

        }

        // Progress pointers.
        linslopeptr++;
        scaleptr += ckd->dim_nonlin_exptime;
        resptr += ckd->dim_nonlin_exptime*ckd->dim_nonlin_spline;
        chi2ptr += ckd->dim_nonlin_exptime;

    }
    int percentage_end = 100*((iexptime+1)*ipix_end+(ckd->dim_nonlin_exptime-iexptime-1)*ipix_start)/(ckd->npix*ckd->dim_nonlin_exptime);
    percentagelog_close(percentage_end);

    return 0;

} // }}}

// Interpolate to chosen set to use.
int Nonlincal::process_finalize( // {{{
)
{

    return 0;

} // }}}

int Nonlincal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{
    for (size_t iexptime=0 ; iexptime<ckd->dim_nonlin_exptime ; iexptime++) {
        DetailedOutputNonlin &det_out = detailed_output_series[iexptime];
        string exptime_suffix = format("_%zu",iexptime);
        string nm = format("nmeas%s",exptime_suffix.c_str());
        NcDim dimid_nmeas;
        netcdf_check(nc,dimid_nmeas = grp.addDim(nm.c_str(),det_out.nl1a));
        vector<NcDim> dims = {ckd->dimid_detector_spat,ckd->dimid_detector_spec,dimid_nmeas};
        nm = format("signal_meas%s",exptime_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dims).putVar(det_out.signal_meas.data()));
        nm = format("signal_corr%s",exptime_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dims).putVar(det_out.signal_corr.data()));
    }

    return 0;

} // }}}

