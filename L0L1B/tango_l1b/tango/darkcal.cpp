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
#include "processor.h"
#include "darkcal.h"

// Settings functions.
Settings_darkcal::Settings_darkcal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "dark";
    output_level_max = 1;
} // }}}
Settings_darkcal::~Settings_darkcal() {} // Destructor.
int Settings_darkcal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(order); // B-spline order for temperature fit (one higher than polynomial order).
    recognize_setting(nominal_temperature); // Standard temperature.
    recognize_setting(saturation_fraction_start); // Fraction of maximum signal to include in initial fit.
    recognize_setting(chi2_increment_threshold); // Threshold of chi square increment where to stop.
    recognize_setting(mask_chi2_max); // Maximum fit residual chi squared for dark fit.
    recognize_setting(mask_offset_min); // Minimum fixed dark signal (independent of integration time).
    recognize_setting(mask_offset_max); // Maximum fixed dark signal (independent of integration time).
    recognize_setting(mask_current_min); // Minimum dark signal added per second of integration time.
    recognize_setting(mask_current_max); // Maximum dark signal added per second of integration time.

    return 0;

} // }}}

// Constructor for the Darkcal CKD structure.
Darkcal::Darkcal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("dark");
    set = make_unique<Settings_darkcal>(this);
    Processor::set = set.get();
    own_skip = &ckd->dark_skip;
} // }}}
Darkcal::~Darkcal() {}

// Check the settings before reading any L1A file.
int Darkcal::process_init( // {{{
)
{

    check_error(set->order < 1,"Error: B-spline degrees of freedom for temperature must be at least one, for a constant over temperature. Setting 'order'");
    check_error(set->nominal_temperature == NC_FILL_DOUBLE,"Error: Missing nominal temperature. Setting 'nominal_temperature'");
    check_error(set->saturation_fraction_start == NC_FILL_DOUBLE,"Error: Missing starting saturation fraction for saturation protection. To use all measurements, use a value higher than one. Setting 'saturation_fraction_start'");
    check_error(set->saturation_fraction_start <= 1.0 && set->chi2_increment_threshold == NC_FILL_DOUBLE,"Error: Missing chi square increment threshold for saturatino protection. Use zero to just use the selected saturation fraction. Setting 'chi2_increment_threshold'");

    return 0;

} // }}}

// Dark correction calibration software.
// This is the planned protocol.
// 1. Extract necessary metadata from input.
// 2. Fit temperature-dependent dark signal.
// 3. Add a noise estimate.
int Darkcal::process_batch( // {{{
    size_t ibatch // Meaningless zero here.
)
{

    vector<double> time(nl1a); // Should end up with size ntime.
    vector<double> temp(nl1a); // Should end up with size ntemp.

    for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
        L1A *l1a = l1a_instances[il1a];
        time[il1a] = l1a->exposure_time;
    }

    // Everything is independent of pixels what so ever.
    vector<double> knots_temp = {set->nominal_temperature,set->nominal_temperature+1}; // Same idea, such that first term is at nominal temperature.
    Bspline b_temp(this,set->order,2,knots_temp.data());
    ckd->dim_dark_order = b_temp.nspline;
    size_t nstate = 2*b_temp.nspline;

    // Construct full Jacobian. The state will be just the temperature B-splines followed
    // by the temperature B-splines multiplied with the exposure time. So the
    // first nspline components are for the dark offset and the last nspline
    // components are for the dark current.
    vector<double> jac(nstate*nl1a);
    handle(b_temp.jaccalc(nl1a,temp.data(),jac.data()));
    // That is the first half of the Jacobian. Now, the second.
    for (size_t ispline=0 ; ispline<b_temp.nspline ; ispline++) {
        for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
            jac[(ispline+b_temp.nspline)*nl1a+il1a] = jac[ispline*nl1a+il1a] * time[il1a];
        }
    }

    // Saturation protection: Figure out what percentage of well filling each
    // measurement is. That will be based on the medians of the full images.
    // Next, a fit is performed as normal, using these medians and based on the
    // chi squared, it is figured out which measurements to include.
    // We could use the same measurements for each pixel, or figure out the
    // subset again for each pixel.
    medians.resize(nl1a); // This is a member variable, because if detailed output is selected, this is included in the output.
    percentagelog_open("Acquire image medians");
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
        percentagelog_progress(il1a,nl1a);
        median(ckd->npix,l1a_instances[il1a]->image,medians[il1a]);
    }
    percentagelog_close();
    // Sort this thing.
    vector<size_t> indexarray(nl1a);
    vector<double> medians_sorted(nl1a);
    Array<double>::sort(nl1a,medians.data(),medians_sorted.data(),indexarray.data());

    writelog(log_trace,"Construct saturation filter based on fit through medians.");
    size_t nl1a_include_start = 0;
    // Include all measurements with less well filling than threshold from setting.
    double mn = medians_sorted[0];
    double mx = medians_sorted[nl1a-1];
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
        double &med = medians_sorted[il1a];
        if ((med-mn) / (mx-mn) < set->saturation_fraction_start) nl1a_include_start++;
    }

    // The fit will fail if there are not enough included L1A images, but
    // if there are just enough points, the fit will be perfect and the
    // chi squared will be undefined, so the measurement selection protocol
    // will fail.
    check_error(nl1a_include_start <= nstate,"Insufficient measurements in initial reference fraction for protocol. Increase setting 'saturation_fraction_start' or include more low measurements.");

    // Create mask.
    vector<bool> mask(nl1a,false);

    // Remove excluded measurements from the Jacobian.
    for (size_t irem=nl1a_include_start ; irem<nl1a ; irem++) {
        mask[indexarray[irem]] = true;
    }
    // Perform a reference retrieval for the chosen start threshold.
    vector<double> res(nstate);
    check_error(linear_invert(nstate,nl1a,jac.data(),OPT_NONE,NULL,&mask,medians.data(),res.data()) != 0,"Error: Linear inversion of medians failed at starting threshold.");
    vector<double> mod(nl1a);
    Matrix::matmul_rl_fold_slow(nstate,nl1a,jac.data(),res.data(),mod.data());
    double chi2_ref = 0.0;
    for (size_t iinc=0 ; iinc<nl1a_include_start ; iinc++) {
        size_t &idx = indexarray[iinc];
        chi2_ref += pow(medians[idx]-mod[idx],2.0) / (nl1a_include_start-nstate);
    }

    if (set->output_level >= 1) {
        det_nl1a = nl1a;
        // Medians need no initialization, because it is also used as working variable.
        det_medianfit_chi2.resize(nl1a,NC_FILL_DOUBLE);
        det_inclusion.resize(nl1a,0); // Defaults to zero (not even considered).
        det_medianfit_chi2[indexarray[nl1a_include_start-1]] = chi2_ref;
        // Mark ones for the set that is included from the start.
        for (size_t il1a=0 ; il1a<nl1a_include_start ; il1a++) det_inclusion[indexarray[il1a]] = 1;
    }

    // If the loop is completed without break, all measurements should be
    // included. Then, nl1a_include will not be overwritten.
    size_t nl1a_include = nl1a;
    for (size_t iinc_add=nl1a_include_start ; iinc_add<nl1a ; iinc_add++) {
        // Restore measurement number iinc_add.
        mask[indexarray[iinc_add]] = false;
        // Perform the retrieval.
        vector<double> res(nstate);
        check_error(linear_invert(nstate,nl1a,jac.data(),OPT_NONE,NULL,&mask,medians.data(),res.data()) != 0,"Error in linear inversion when relaxing threshold with medians.");
        vector<double> mod(nl1a);
        Matrix::matmul_rl_fold_slow(nstate,nl1a,jac.data(),res.data(),mod.data());
        double chi2_test = 0.0;
        // Index iinc_add is included now.
        for (size_t iinc=0 ; iinc<=iinc_add ; iinc++) {
            size_t &idx = indexarray[iinc];
            chi2_test += pow(medians[idx]-mod[idx],2.0) / (nl1a_include_start-nstate);
        }
        if (set->output_level >= 1) {
            det_medianfit_chi2[indexarray[iinc_add]] = chi2_test;
            det_inclusion[indexarray[iinc_add]] = iinc_add - nl1a_include_start + 2; // If filter does not trigger. Will be overwritten if filter triggers.
        }
        if (chi2_test > set->chi2_increment_threshold*chi2_ref) {
            // Remove last measurement again.
            mask[indexarray[iinc_add]] = true;
            // Now, index iinc_add is not included anymore, so total number is equal to iinc_add.
            nl1a_include = iinc_add; // Overwrite this thing.
            if (set->output_level >= 1) det_inclusion[indexarray[iinc_add]] = -1; // Triggered filter.
            break;
        }
    }
    writelog(log_debug,"Including %zu measurements in the fit.",nl1a_include);
    // Regardless of whether the loop is completed or broken. nl1a_include is
    // now set to the right value.
    // More importantly, mask is now correct.

    // We make a gain matrix that is complete, it will have zeros for the
    // excluded measurements. That has to do with the misery that would arise
    // because they are not sorted.
    vector<double> gain(nstate*nl1a);
    check_error(linear_invert_gain(nstate,nl1a,jac.data(),OPT_NONE,NULL,&mask,gain.data()) != 0,"Error while calculating dark correction gain matrix.");

    // Construct output.
    ckd->dark_offset.resize(ckd->dim_dark_order*ckd->npix);
    ckd->dark_current.resize(ckd->dim_dark_order*ckd->npix);
    dark_chi2.resize(ckd->npix);

    // Moving pointers.
    double *dark_offset_cur = ckd->dark_offset.data();
    double *dark_current_cur = ckd->dark_current.data();
    double *dark_chi2_cur = dark_chi2.data();

    percentagelog_open("Fitting temperature-dependent dark signal per pixel");
    for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {

        percentagelog_progress(ipix,ckd->npix);

        vector<double> meas(nl1a);
        for (size_t il1a=0 ; il1a<nl1a ; il1a++) meas[il1a] = l1a_instances[il1a]->image[ipix];
        vector<double> res(nstate);

        Matrix::matmul_rl_fold_quick(nstate,nl1a,gain.data(),meas.data(),res.data());

        // Split results into dark offset and dark current.
        memcpy(dark_offset_cur,res.data(),b_temp.nspline*sizeof(double));
        memcpy(dark_current_cur,&res[b_temp.nspline],b_temp.nspline*sizeof(double));

        // Calculate chi squared.
        *dark_chi2_cur = 0.0;
        vector<double> mod(nl1a);
        Matrix::matmul_rl_fold_slow(nstate,nl1a,jac.data(),res.data(),mod.data());
        for (size_t iinc=0 ; iinc<nl1a_include ; iinc++) {
            size_t &idx = indexarray[iinc];
            *dark_chi2_cur += pow(meas[idx]-mod[idx],2.0) / (nl1a_include - nstate);
        }

        // Apply pixel mask.
        if (set->mask_chi2_max != NC_FILL_DOUBLE && *dark_chi2_cur > set->mask_chi2_max) ckd->mask[ipix] = true;
        // The offset and current is checked for nominal temperature.
        // That is the first element.
        if (set->mask_offset_min != NC_FILL_DOUBLE && dark_offset_cur[0] < set->mask_offset_min) ckd->mask[ipix] = true;
        if (set->mask_offset_max != NC_FILL_DOUBLE && dark_offset_cur[0] > set->mask_offset_max) ckd->mask[ipix] = true;
        if (set->mask_current_min != NC_FILL_DOUBLE && dark_current_cur[0] < set->mask_current_min) ckd->mask[ipix] = true;
        if (set->mask_current_max != NC_FILL_DOUBLE && dark_current_cur[0] > set->mask_current_max) ckd->mask[ipix] = true;

        // Progress pointers.
        dark_offset_cur += ckd->dim_dark_order;
        dark_current_cur += ckd->dim_dark_order;
        dark_chi2_cur++;

    }
    percentagelog_close();

    // Do not forget to communicate the nominal temperature to the CKD.
    ckd->dark_nominal_temperature = set->nominal_temperature;

    // }}}

    // Set diagnostic CKD pointer.
    ckd->diag_dark_chi2 = dark_chi2.data();

    return 0;

} // }}}

int Darkcal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{

    // Dimension should still exist as nl1a and nl1a_total. But if anything
    // is changed in batching, this becomes unreliable, so used the saved
    // variable.
    NcDim dimid_nmeas;
    netcdf_check(nc,dimid_nmeas = grp.addDim("nmeas",det_nl1a));
    netcdf_check(nc,grp.addVar("medians",ncDouble,dimid_nmeas).putVar(medians.data()));
    netcdf_check(nc,grp.addVar("medianfit_chi2",ncDouble,dimid_nmeas).putVar(det_medianfit_chi2.data()));
    NcVar var_inclusion;
    netcdf_check(nc,var_inclusion = grp.addVar("inclusion",ncInt,dimid_nmeas));
    netcdf_check(nc,var_inclusion.putAtt("meaning","1: In initial set. >1: Iteratively included and survived filter (in order 2,n). -1: Triggered filter, thus excluded. 0: Higher than measurement that triggered filter, so not even considered."));
    netcdf_check(nc,var_inclusion.putVar(det_inclusion.data()));

    return 0;

} // }}}

