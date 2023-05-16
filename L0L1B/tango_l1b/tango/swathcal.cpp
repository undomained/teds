// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "array.h"
#include "matrix.h"
#include "lininv.h"
#include "vector.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "splitgaussfit.h"
#include "ckd.h"
#include "l1a.h"
#include "theodolite.h"
#include "batch.h"
#include "processor.h"
#include "swathcal.h"

// Settings functions.
Settings_swathcal::Settings_swathcal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "swath";
    output_level_max = 1;
    opt.nonlin = false;
    opt.stray = false;
} // }}}
Settings_swathcal::~Settings_swathcal() {}
int Settings_swathcal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(theodolite_reference_filename); // NetCDF file name of the theodolite reference file (a lot of measurement on the wrong day).
    recognize_setting(theodolite_spotcheck_filename); // NetCDF file name of the theodolite spot-check file (one measurement on correct day).
    recognize_setting_vector(viewport_alt_angles); // Expected ALT angles of each viewport (also used to put the viewports in the correct order).
    recognize_setting(theodolite_alt_angle_tol); // Maximum allowed different between ALT angle of theodolite measurement and expected ALT angle of a viewport.
    recognize_setting(act_linesearch_tol); // Convergence tolerance in line search over ACT angles.
    recognize_setting(act_angle_tol); // Difference between across-track angles that are considered to be the same.
    recognize_setting(actscan_sampling); // Sampling for initial scan over swath to find the correct FOVcal ACT angle (default = just use swathcal stage ACT angle).
    recognize_setting(insignificant_island); // Setting for peak detection: Size of peaks or non-peaks to be discarded.
    recognize_setting(maxreso); // Highest denominator in setting threshold to detect peaks.
    recognize_setting(backgroundorder); // Degrees of freedom of background polynomial for peak detection.
    recognize_setting(gaussfit); // Flag for performing gauss fit.
    recognize_setting(domain_range); // Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    recognize_setting(maxiter); // Number of iterations for fitting Gauss peak shape.
    recognize_setting(magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
    recognize_setting(postol); // Tolerance of peak positions (in indices).
    recognize_setting(widthtol); // Tolerance of peak FWHM (in indices).
    recognize_setting(backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
    recognize_setting(steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
    recognize_setting(initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
    recognize_setting(reducerfactor_fail); // Multiplication of step reducer for a failed step.
    recognize_setting(reducerfactor_success); // Multiplication of step reducer for a successful step.
    recognize_setting(reducer_limit); // Maximum step control reducer for convergence.
    recognize_setting_vector(order); // Degrees of freedom for polynomial fit of ALT against ACT angles per viewport, each with a minimum number of successful associations.
    recognize_setting_vector(min_fit); // Minimum required associations to perform the fit with the respective order.
    recognize_setting_vector(rotmat); // Rotation matrix from spacecraft to instrument coordinates. Nine elements. Source direction is quick. Destination dimension is slow.
    recognize_setting(rot_warningthreshold); // Allowed error in check if rotmat is a proper rotation matrix.
    recognize_setting(theo_altmodel_maxiter); // Maximum number of iterations in non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_tol); // Convergence tolerance on vector elements in non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_steptolerance); // Inversion step acceptance criterion for non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter for non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_reducerfactor_fail); // Multiplication of step reducer for a failed step for non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_reducerfactor_success); // Multiplication of step reducer for a successful step for non-linear inversion for ALT-model in the theodolite.
    recognize_setting(theo_altmodel_reducer_limit); // Maximum step control reducer for convergence for non-linear inversion for ALT-model in the theodolite.

    return 0;

} // }}}

// Constructor for the Swathcal CKD structure.
Swathcal::Swathcal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("swath");
    set = make_unique<Settings_swathcal>(this);
    Processor::set = set.get();
    own_skip = &ckd->swath_skip;
} // }}}
Swathcal::~Swathcal() {}

// Swath calculation.
// This is a placeholder. It should calculate the swath vectors from calibration
// measurements, but now, it just makes up swath vectors.
int Swathcal::process_init( // {{{
)
{

    // Check legal settings.
    check_error(set->viewport_alt_angles.size() != ckd->dim_vp,"Error: Wrong number (%zu) of expected ALT angles given. One per viewport (%zu) should be given. Setting 'viewport_alt_angles'",set->viewport_alt_angles.size(),ckd->dim_vp);
    check_error(set->theodolite_alt_angle_tol < 0,"Error: Missing or negative tolerance for along-track stage angle for sorting theodolite measurements. Setting 'theodolite_alt_angle_tol'");
    // Check if any tolerance domain overlaps between different viewports.
    for (size_t ivp1=0 ; ivp1<ckd->dim_vp-1 ; ivp1++) {
        for (size_t ivp2=ivp1+1 ; ivp2<ckd->dim_vp ; ivp2++) {
            double diff = abs(set->viewport_alt_angles[ivp2]-set->viewport_alt_angles[ivp1]);
            check_error(diff < 2.0*set->theodolite_alt_angle_tol,"Error: With these expected ALT angles and tolerance value, viewports %zu and %zu overlap.",ivp1,ivp2);
        }
    }
    check_error(set->act_angle_tol < 0,"Error: Missing or negative across-track angle tolerance. Setting 'act_angle_tol'");
    check_error(set->act_linesearch_tol <= 0.0,"Error: Zero, negative or missing convergence tolerance for ACT angle line search. Setting 'act_linesearch_tol'");
    check_error(set->actscan_sampling != NC_FILL_DOUBLE && set->actscan_sampling <= 0.0,"Error: Zero or negative sampling for pre-scan over ACT angles for detector positioning algorithm. Setting 'actscan_sampling'");
    check_error(set->maxreso < 2,"Error: Maximum peak-detection threshold resolution must be at least two. The higher this value (only powers of two make sense), the more patience the program has to find peaks. Setting 'maxreso'");
    check_error(set->backgroundorder < 1,"Error: B-spline degrees of freedom for the background for the Gauss fit must be at least one, for a constant background. Setting 'backgroundorder'");
    // More inversion settings. Only relevant for a gauss fit.
    if (set->gaussfit) {
        check_error(set->domain_range == NC_FILL_DOUBLE,"Error: Missing domain range to determine fitting domain for the Guass, necessary because gaussfit is turned on. Setting 'domain_range'");
        check_error(set->maxiter < 1,"Error: Missing or zero maximum iterations in non-linear inversion. Setting 'maxiter'");
        check_error(set->magtol == NC_FILL_DOUBLE,"Error: Missing Gauss magnitude tolerance for convergence. Setting 'magtol'");
        check_error(set->postol == NC_FILL_DOUBLE,"Error: Missing Gauss position tolerance for convergence. Setting 'postol'");
        check_error(set->widthtol == NC_FILL_DOUBLE,"Error: Missing Gauss FWHM tolerance for convergence. Setting 'widthtol'");
        check_error(set->backgroundtol == NC_FILL_DOUBLE,"Error: Missing background tolerance in Gauss fit. Setting 'backgroundtol'");
        check_error(set->steptolerance < 1.0,"Error: Inversion step tolerance below one may lead to an endless loop. Setting 'steptolerance'");
        check_error(set->initial_step_reducer == NC_FILL_DOUBLE,"Error: Missing initial Levenber-Marquardt step reducer. Setting 'initial_step_reducer'");
        if (set->initial_step_reducer != 0.0) {
            check_error(set->reducerfactor_success == NC_FILL_DOUBLE,"Error: Missing step reducer factor for successful steps. Settings 'reducerfactor_success'");
            check_error(set->reducerfactor_success > 1.0,"Error: Step reduction step for successful steps should be not higher than one. Settings 'reducerfactor_success'");
        }
        if (set->steptolerance == NC_FILL_DOUBLE) {
            // Step failure is impossible.
            if (set->reducerfactor_fail != NC_FILL_DOUBLE) writelog(log_warning,"Warning: Step reducer factor for failed steps ignored, because no step tolerance is given, so steps cannot be rejected.");
        } else {
            // Check legal step reduction.
            check_error(set->reducerfactor_fail == NC_FILL_DOUBLE,"Error: Missing reduction factor for failed steps. Stting 'reducerfactor_fail'");
            check_error(set->reducerfactor_fail <= 1.0,"Error: Step reduction factor for failed steps should be higher than one. Stting 'reducerfactor_fail'");
        }
    }
    nattempt = set->min_fit.size();
    size_t dimcheck = set->order.size();
    check_error(dimcheck != nattempt,"Settings error: Array lengths of 'order' and 'min_fit' should be equal. Now, 'order' is %zu large and 'min_fit' is %zu large.",dimcheck,nattempt);
    for (size_t iattempt=0 ; iattempt<nattempt ; iattempt++) {
        check_error(set->order[iattempt] > set->min_fit[iattempt],"Error: Minimum number of fits (%zu) too low for a fit with DFS (%zu).",set->min_fit[iattempt],set->order[iattempt]);
    }
    for (size_t iattempt=1 ; iattempt<nattempt ; iattempt++) {
        check_error(set->min_fit[iattempt] >= set->min_fit[iattempt-1],"Error: Minimum required successful fits should decrease per attempt. Setting 'min_fit'");
    }
    // Rotation matrix.
    check_error(set->rotmat.size() != 9,"Error: Rotation matrix should have nine elements. Setting 'rotmat'");
    vector<double> rtr(6);
    Matrix::matmul_rr_fold_quick_sym(3,3,set->rotmat.data(),set->rotmat.data(),rtr.data());
    if (abs(rtr[0]-1.0) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of one at location (X,X)",rtr[0]);
    if (abs(rtr[1]) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of zero at location (X,Y) and (Y,X)",rtr[1]);
    if (abs(rtr[2]-1.0) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of one at location (Y,Y)",rtr[2]);
    if (abs(rtr[3]) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of zero at location (X,Z) and (Z,X)",rtr[3]);
    if (abs(rtr[4]) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of zero at location (Y,Z) and (Z,Y)",rtr[4]);
    if (abs(rtr[5]-1.0) > set->rot_warningthreshold) writelog(log_warning,"Warning: Rotation matrix check times its transpose has %.12f instead of one at location (Z,Z)",rtr[5]);
    double det = 
        set->rotmat[0]*set->rotmat[4]*set->rotmat[8] +
        set->rotmat[1]*set->rotmat[5]*set->rotmat[6] +
        set->rotmat[2]*set->rotmat[3]*set->rotmat[7] -
        set->rotmat[0]*set->rotmat[5]*set->rotmat[7] -
        set->rotmat[1]*set->rotmat[3]*set->rotmat[8] - 
        set->rotmat[2]*set->rotmat[4]*set->rotmat[6]
    ;
    check_error(det < 0.0,"Error: Rotation matrix has negative (%.7f) determinant. Possibly, it is a roto-reflection. Setting 'rotmat'.",det);

    // Transform to quaternion so that rotation can be applied in one line.
    double realpart = 0.5*sqrt(set->rotmat[0]+set->rotmat[4]+set->rotmat[8]+1.0);
    q_inst = Quaternion(
        (set->rotmat[7]-set->rotmat[5]) / (4.0*realpart),
        (set->rotmat[2]-set->rotmat[6]) / (4.0*realpart),
        (set->rotmat[3]-set->rotmat[1]) / (4.0*realpart),
        realpart
    ).conjug();
    // Input is conjugated. The input is the spacecraft to SPEX matrix.
    // We use a quaternion to rotate SPEX to the spacecraft.
    q_inst.normalize(); // Getting rid of eventual rounding errors possibly within the warning threshold.

    // Log the instrument rotation.
    if (q_inst[3] == 1.0) {
        writelog(log_debug,"No instrument rotation.");
    } else {
        double ang = 2.0*acos(q_inst[3]);
        Vector axis(q_inst[0],q_inst[1],q_inst[2]);
        axis.normalize();
        writelog(log_debug,"Instrument rotation %.7f degrees along axis (%.7f,%.7f,%.7f)",ang/DEGREES,axis[0],axis[1],axis[2]);
    }

    // Convert units of settings.
    set->act_angle_tol *= DEGREES;
    set->act_linesearch_tol *= DEGREES;
    if (set->actscan_sampling != NC_FILL_DOUBLE) set->actscan_sampling *= DEGREES;
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        set->viewport_alt_angles[ivp] *= DEGREES;
    }
    set->theodolite_alt_angle_tol *= DEGREES;

    // Make batches per viewport per ACT angle.

    // The L1A format allows illumination of multiple viewport, but these
    // should all have a star stimulus, so that should not be possible.
    // Therefore, we assume and assert that all files only have one
    // illuminated viepowrt (and also not zero).

    // Minimum number of ACT angles per viewport to make a fit.
    min_fit_poorest = set->min_fit[nattempt-1];

    // Allocate the batches array with more then enough batches.
    batches.resize(nl1a_total);
    nbatch = 0; // Counter.
    // Running batch pointer.
    Batch *bat_cur = batches.data();

    // Administration.
    ibatch_start.resize(ckd->dim_vp);
    nbatch_vp.resize(ckd->dim_vp);
    nfail_vp.resize(ckd->dim_vp,0);

    // Sort the ACT angles.
    vector<double> act_angles_unsorted(nl1a_total);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) act_angles_unsorted[il1a] = l1a_instances_all[il1a]->act_angle;
    vector<double> act_angles_sorted(nl1a_total);
    vector<size_t> l1a_indexarray(nl1a_total);
    Array<double>::sort(nl1a_total,act_angles_unsorted.data(),act_angles_sorted.data(),l1a_indexarray.data());

    // I do not want to use pow(2,ivp), because pow is a float function.
    // So, I will iterate.
    uint8_t l1a_viewport_cur = 1;
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        ibatch_start[ivp] = nbatch;
        if (ckd->vp_mask[ivp]) {
            nbatch_vp[ivp] = 0;
            l1a_viewport_cur *= 2;
            continue;
        }
        double *act_prv = NULL; // Previous ACT angle.
        // Open first empty batch in the viewport.
        bat_cur->setSize(nl1a_total); // Too large.
        size_t nl1a_batch = 0;
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            L1A *l1a = l1a_instances_all[l1a_indexarray[il1a]].get();
            if (l1a->viewport != l1a_viewport_cur) continue; // Only consider own viewport.
            double *act_cur = &act_angles_sorted[il1a];
            if (act_prv == NULL || (*act_cur - *act_prv <= set->act_angle_tol)) {
                // Write L1A instance in current batch.
                bat_cur->l1a[nl1a_batch] = l1a;
                nl1a_batch++;
            } else {
                // Close off previous batch.
                bat_cur->setSize(nl1a_batch); // Correct size.
                bat_cur++;
                nbatch++;
                // Open new batch with this L1A in it.
                bat_cur->setSize(nl1a_total); // Too large.
                nl1a_batch = 0;
                // Write file in new batch.
                bat_cur->l1a[nl1a_batch] = l1a;
                nl1a_batch++;
            }
            act_prv = act_cur;
        }
        // Close off last batch in the viewport.
        check_error(nl1a_batch == 0,"Error: There seems to be no L1A file for viewport %zu.",ivp);
        bat_cur->setSize(nl1a_batch);
        bat_cur++;
        nbatch++;
        // Administration.
        nbatch_vp[ivp] = nbatch-ibatch_start[ivp];
        check_error(nbatch_vp[ivp] < min_fit_poorest,"Error: Viewport %zu only has %zu different ACT angles, while the simplest fit already requires %zu successful fits.",ivp,nbatch_vp[ivp],min_fit_poorest);
        l1a_viewport_cur *= 2;
    }

    // Encapsulate batch array.
    batches.resize(nbatch);

    // Check that all L1A files are from one viewport only.
    // If there are skipped viewports, not all L1A files are in a batch. That is okay.
    size_t cnt = 0;
    l1a_viewport_cur = 1;
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            L1A *l1a = l1a_instances_all[il1a].get();
            if (l1a->viewport == l1a_viewport_cur) cnt++;
        }
        l1a_viewport_cur *= 2;
    }
    check_error(cnt != nl1a_total,"Error: Not all L1A files are from exactly one viewport.");

    // Shape arrays.
    act_fovcal.resize(nbatch); // Used to express polynomials in and sample at L1B FOVs (based on detector positions).
    act_swathcal.resize(nbatch); // Stage ACT angles during these measurements.
    alt_swathcal.resize(nbatch); // Peak ALT angles.
    // Both act_swathcal and alt_swathcal are expressed as polynomials of act_fovcal.
    // Here, act_fovcal represents the detector positioning.
    // They are sampled on the act_fovcals from the L1B (are in CKD).
    // And then, the swathcal stage angles are fed to the theodolite
    // system and turned into swath vectors.
    mask.resize(nbatch); // Mask, for if an ACT angle seems not to.
    mag.resize(nbatch); // Magnitude of fit, for weighing or filtering.
    // During finalize, the fit will be performed.

    theo = make_unique<Theodolite>(this);
    handle(theo->calculate_reference(
        ckd->dim_vp,
        set->theodolite_reference_filename,
        set->viewport_alt_angles.data(),
        set->theodolite_alt_angle_tol,
        set->theo_altmodel_maxiter,
        set->theo_altmodel_tol,
        set->theo_altmodel_steptolerance,
        set->theo_altmodel_initial_step_reducer,
        set->theo_altmodel_reducerfactor_fail,
        set->theo_altmodel_reducerfactor_success,
        set->theo_altmodel_reducer_limit
    ));

    handle(theo->calculate_spotcheck(
        set->theodolite_spotcheck_filename
    ));

    return 0;

} // }}}

int Swathcal::process_batch( // {{{
    size_t ibatch // This is the viewport index.
)
{


    // Figure out which viewport it is.
    size_t ivp = ckd->dim_vp; // Will be overwritten.
    for (size_t attempt=0 ; attempt<ckd->dim_vp ; attempt++) {
        if (ibatch >= ibatch_start[attempt]) ivp = attempt;
        // Do not break. The last viewport in which this condition is true
        // is the right one.
    }
    check_error(ivp == ckd->dim_vp,"Program error: This should never happen. Debug code using traceback.");
    check_error(ckd->vp_mask[ivp],"Program error: Batch associated with skipped viewport. That should never be possible.");

    // If anything in the process fails, we will throw away the entire
    // batch. That means that eventual results in the middle will also
    // be discarded.
    bool fail = true;
    // Breakable do-once loop.
    for (bool ctd=true ; ctd ; ctd=false) {

        // Write swathcal ACT angles. They are all the same inside the batch
        // besides eventual rounding differences covered by a tolerance.
        act_swathcal[ibatch] = 0.5*(l1a_instances[0]->act_angle + l1a_instances[nl1a-1]->act_angle); // In principle, they are all equal.
        writelog(log_info,"Viewport %zu, ACT angle %.7f",ivp,act_swathcal[ibatch]/DEGREES);

        // Figure out where the spectrum is on the detector. By doing so,
        // the corresponding FOVcal ACT angle is derived, which because
        // of renewal of the laboratory setup may differ from the ACT angle
        // at the measruements themselves (swathcal). To ensure that there
        // is a spectrum, we add up all the images as they should have their
        // either-or-not visible spectrum on the same detector pixels.
        vector<double> image_agg(ckd->npix,0.0);
        vector<bool> pixelmask_agg(ckd->npix,false);

        // Add L1As to the image and mask, but ignore the noise.
        for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
            // This only works for non-binned L1A. It can be extended so
            // that it works with L1As with all the same binning table.
            for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
                if (l1a_instances[il1a]->pixelmask[ipix]) pixelmask_agg[ipix] = true;
                image_agg[ipix] += l1a_instances[il1a]->image[ipix];
            }
        }

        vector<double> linesearch_act(3,NC_FILL_DOUBLE);
        double linesearch_strength = NC_FILL_DOUBLE;
        double &act_swath_start = ckd->fov_act_angles[0];
        double &act_swath_end = ckd->fov_act_angles[ckd->fov_nfov_vp[ivp]-1];
        if (set->actscan_sampling == NC_FILL_DOUBLE) {
            // If the swathcal stage ACT angle is outside the domain
            // of the swath, this is a wrong a-priori and the fit cannot
            // be done.
            if (act_swathcal[ibatch] < act_swath_start || act_swathcal[ibatch] > act_swath_end) break; // Fail. The a-priori is not in the swath range.
            linesearch_act[0] = act_swath_start;
            linesearch_act[1] = act_swathcal[ibatch];
            linesearch_act[2] = act_swath_end;
            double linesearch_strength_left = NC_FILL_DOUBLE;
            double linesearch_strength_right = NC_FILL_DOUBLE;
            handle(act_extract(ivp,image_agg.data(),pixelmask_agg,linesearch_act[0],linesearch_strength_left)); // Last argument is output.
            handle(act_extract(ivp,image_agg.data(),pixelmask_agg,linesearch_act[1],linesearch_strength)); // Last argument is output.
            handle(act_extract(ivp,image_agg.data(),pixelmask_agg,linesearch_act[2],linesearch_strength_right)); // Last argument is output.
            // Only accept this if the a-priori is the best of the three.
            if (linesearch_strength_left > linesearch_strength || linesearch_strength_right > linesearch_strength) break; // Fail. One of the edges of the swath are stronger than the a-priori.
        } else {
            // Slightly change this sampling so that an integer number of
            // samples fit in the swath.
            double n_interval = (act_swath_end-act_swath_start) / set->actscan_sampling;
            size_t nsamp = static_cast<size_t>(n_interval) + 2; // Plus one for taking the ceil instead of the floor, plus one to convert interval to samples.
            check_error(nsamp == 2,"Error: Too coarse scan sampling. No intermediate angles are sampled.");
            vector<double> act_scan(nsamp,NC_FILL_DOUBLE);
            vector<double> strength_scan(nsamp,NC_FILL_DOUBLE);
            for (size_t isamp=0 ; isamp<nsamp ; isamp++) {
                act_scan[isamp] = (isamp*act_swath_end + (nsamp-isamp-1)*act_swath_start) / (nsamp-1);
                handle(act_extract(ivp,image_agg.data(),pixelmask_agg,act_scan[isamp],strength_scan[isamp])); // Last argument is output.
            }
            // Find the peak.
            double mx = strength_scan[0];
            size_t idx_mx = 0;
            for (size_t isamp=1 ; isamp<nsamp ; isamp++) {
                double &newstrength = strength_scan[isamp];
                if (newstrength > mx) {
                    mx = newstrength;
                    idx_mx = isamp;
                }
            }
            if (idx_mx == 0 || idx_mx == nsamp-1) {
                // Exception, peak way on the left or right.
                bool left = idx_mx == 0;
                size_t idx_side;
                size_t idx_next_to_side;
                size_t idx_ls_side;
                size_t idx_ls_space;
                if (left) {
                    idx_side = 0;
                    idx_next_to_side = 1;
                    idx_ls_side = 0;
                    idx_ls_space = 2;
                } else {
                    idx_side = nsamp-1;
                    idx_next_to_side = nsamp-2;
                    idx_ls_side = 2;
                    idx_ls_space = 0;
                }
                    
                double prv_act = act_scan[idx_next_to_side];
                double new_act = NC_FILL_DOUBLE;
                double new_strength = NC_FILL_DOUBLE;
                while (prv_act - act_scan[idx_side] > set->act_linesearch_tol) {
                    new_act = 0.5*(act_scan[idx_side]+prv_act);
                    handle(act_extract(ivp,image_agg.data(),pixelmask_agg,new_act,new_strength)); // Last argument is output.
                    if (new_strength > strength_scan[idx_side]) break; // Breaks the while-loop. This is not a failure. On the contrary, this will induce success.
                    prv_act = new_act;
                }
                // If the tolerance is reached and no new maximum is found, it is a failure.
                if (new_strength <= strength_scan[idx_side]) break; // Fail.
                linesearch_act[idx_ls_side] = act_scan[idx_side];
                linesearch_act[1] = new_act;
                linesearch_act[idx_ls_space] = prv_act;
                linesearch_strength = new_strength;
            } else {
                linesearch_act[0] = act_scan[idx_mx-1];
                linesearch_act[1] = act_scan[idx_mx];
                linesearch_act[2] = act_scan[idx_mx+1];
                linesearch_strength = strength_scan[idx_mx];
            }
        }
        // Converge the line search.
        while (linesearch_act[2] - linesearch_act[0] > set->act_linesearch_tol) {
            double new_act_left = 0.5*(linesearch_act[0] + linesearch_act[1]);
            double new_act_right = 0.5*(linesearch_act[1] + linesearch_act[2]);
            double new_strength_left = NC_FILL_DOUBLE;
            double new_strength_right = NC_FILL_DOUBLE;
            handle(act_extract(ivp,image_agg.data(),pixelmask_agg,new_act_left,new_strength_left)); // Last argument is output.
            handle(act_extract(ivp,image_agg.data(),pixelmask_agg,new_act_right,new_strength_right)); // Last argument is output.
            if (new_strength_left > linesearch_strength && new_strength_right > linesearch_strength) {
                writelog(log_warning,"During ACT line search, both sides seem to improve, domain %.12f degrees",(linesearch_act[2]-linesearch_act[0])/DEGREES);
            }
            if (new_strength_left > linesearch_strength || new_strength_right > linesearch_strength) {
                if (new_strength_left > new_strength_right) {
                    // Left wins.
                    linesearch_strength = new_strength_left;
                    linesearch_act[2] = linesearch_act[1];
                    linesearch_act[1] = new_act_left;
                    // linesearch_act[0] stays the same.
                } else {
                    // Right wins.
                    linesearch_strength = new_strength_right;
                    linesearch_act[0] = linesearch_act[1];
                    linesearch_act[1] = new_act_right;
                    // linesearch_act[2] stays the same.
                }
            } else {
                // Middle wins.
                linesearch_act[0] = new_act_left;
                linesearch_act[2] = new_act_right;
                // linesearch_act[1] and linesearch_strength stay the same.
            }
        }
        act_fovcal[ibatch] = linesearch_act[1];
        writelog(log_debug,"Effective (FOVcal) ACT angle: %.7f",act_fovcal[ibatch]/DEGREES);

        // Now, the contents of the algorithm. Get the along-track angle.

        // Sort the L1A instances by ALT angle.
        vector<double> alt_angles_unsorted(nl1a);
        for (size_t il1a=0 ; il1a<nl1a ; il1a++) alt_angles_unsorted[il1a] = l1a_instances[il1a]->alt_angle;
        vector<double> alt_angles_sorted(nl1a);
        vector<size_t> l1a_indexarray(nl1a);
        Array<double>::sort(nl1a,alt_angles_unsorted.data(),alt_angles_sorted.data(),l1a_indexarray.data());

        // Signal converted to a scalar. This is a measure of how strong the
        // signal with that ALT angle is.
        vector<double> signalmedians(nl1a);

        for (size_t il1a=0 ; il1a<nl1a ; il1a++) {

            handle(act_extract(
                ivp,
                l1a_instances[l1a_indexarray[il1a]]->image,
                l1a_instances[l1a_indexarray[il1a]]->pixelmask,
                act_fovcal[ibatch],
                signalmedians[il1a]
            ));

        }

        // Now we have the median signals as function of the ALT angle.
        // Abscissa: ALT-angle is alt_angles_sorted.
        // Ordinate: Median signals is signalmedians.

        // Fit a Gauss through these.
        // Peak fitter.
        Splitgaussfit spl(this,1,set->backgroundorder,nl1a); // One peak.
        spl.setGaussfit(set->gaussfit); // Flag for performing gauss fit.
        spl.setMaxiter(set->maxiter); // Number of iterations for fitting Gauss peak shape.
        spl.setMagtol(set->magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
        spl.setPostol(set->postol); // Tolerance of peak positions (in indices).
        spl.setWidthtol(set->widthtol); // Tolerance of peak FWHM (in indices).
        spl.setBackgroundtol(set->backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
        spl.setSteptolerance(set->steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
        spl.setInitialStepReducer(set->initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
        spl.setReducerfactorFail(set->reducerfactor_fail); // Multiplication of step reducer for a failed step.
        spl.setReducerfactorSuccess(set->reducerfactor_success); // Multiplication of step reducer for a successful step.
        spl.setReducerLimit(set->reducer_limit); // Maximum step control reducer for convergence.

        // Set abscissa and ordinate.
        spl.setX(alt_angles_sorted.data());
        spl.setMeas(signalmedians.data());
        // Too complicated to take noise into account.
        // Pixel mask should not be there.

        vector<double> res(3);
        if (spl.solve(set->insignificant_island,set->domain_range,set->maxreso,res.data()) != 0) break; // Fail. No proper peak fitted.
        alt_swathcal[ibatch] = res[1];
        mask[ibatch] = false;
        mag[ibatch] = res[0];
        writelog(log_debug,"Associated ALT angle %.7f.",alt_swathcal[ibatch]/DEGREES,ivp);
        fail = false;

    }
    if (fail) {
        act_fovcal[ibatch] = NC_FILL_DOUBLE;
        alt_swathcal[ibatch] = NC_FILL_DOUBLE;
        mask[ibatch] = true;
        mag[ibatch] = NC_FILL_DOUBLE;
        writelog(log_debug,"Fit failed");
        nfail_vp[ivp]++;
        // Check if this viewport lost its chance to have enough successful fits.
        check_error(nbatch_vp[ivp]-nfail_vp[ivp] < min_fit_poorest,"Error: Too many failed fit in viewport %zu. Now, it is impossible to reach %zu successful fits, which is required for the simplest fit.",ivp,min_fit_poorest);
    }

    return 0;

} // }}}

// Private function for spectrum extraction, an often-repeated action
// during process_batch. This extraction does not work for binned images.
int Swathcal::act_extract( // {{{
    size_t ivp, // Viewport index.
    double *image, // Detector image.
    vector<bool> &pixelmask, // Detector pixel mask.
    double act, // ACT-angle in FOVcal perspective.
    double &strength // Output spectrum strength (median of two spectra).
)
{

    // Interpolate over FOV ACT angles.
    // TODO: Maybe we should do this with a spline, also in the FOV.
    double *ckd_act_cur = &ckd->fov_act_angles[0];
    size_t ifov_left = NC_FILL_UINT64;
    if (act <= ckd_act_cur[0]) {
        ifov_left = 0;
    } else if (act >= ckd_act_cur[ckd->fov_nfov_vp[ivp]-1]) {
        ifov_left = ckd->fov_nfov_vp[ivp]-2;
    } else {
        for (size_t ifov_right = 1 ; ifov_right < ckd->fov_nfov_vp[ivp] ; ifov_right++) {
            if (ckd_act_cur[ifov_right] >= act) {
                ifov_left = ifov_right - 1;
                break;
            }
        }
    }
    check_error(ifov_left == NC_FILL_UINT64,"Program error: Variable ifov_left not properly initialization during search loop.");
    // Get the weight factor. This implies linear extrapolation if needed.
    double weightleft = (ckd_act_cur[ifov_left+1]-act) / (ckd_act_cur[ifov_left+1]-ckd_act_cur[ifov_left]);

    // Interpolate the FOV CKD on this ACT angle.
    double *ckd_ispat_left = &ckd->fov_ispat[(ifov_left)*DIM_POL*ckd->dim_detector_spec];
    double *ckd_ispat_right = &ckd->fov_ispat[(ifov_left+1)*DIM_POL*ckd->dim_detector_spec];

    // Extract the spectra to get the median.
    vector<double> spectra(DIM_POL*ckd->dim_detector_spec);
    for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
        for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
            double fov_ispat_interpol = weightleft * ckd_ispat_left[ipol*ckd->dim_detector_spec+ispec] + (1.0-weightleft) * ckd_ispat_right[ipol*ckd->dim_detector_spec+ispec];
            int ispat_left = static_cast<int>(fov_ispat_interpol); // Floor.
            double weightleft = ispat_left+1-fov_ispat_interpol;
            size_t ipix_left = ispat_left*ckd->dim_detector_spec + ispec;
            size_t ipix_right = (ispat_left+1)*ckd->dim_detector_spec + ispec;
            if (ispat_left < 0 || ispat_left > (int)ckd->dim_detector_spat-2) {
                spectra[ipol*ckd->dim_detector_spec+ispec] = 0.0;
            } else {
                if (pixelmask[ipix_left]) {
                    if (pixelmask[ipix_right]) {
                        spectra[ipol*ckd->dim_detector_spec+ispec] = 0.0;
                    } else {
                        spectra[ipol*ckd->dim_detector_spec+ispec] = image[ipix_right];
                    }
                } else {
                    if (pixelmask[ipix_right]) {
                        spectra[ipol*ckd->dim_detector_spec+ispec] = image[ipix_left];
                    } else {
                        spectra[ipol*ckd->dim_detector_spec+ispec] = weightleft * image[ipix_left] + (1.0-weightleft) * image[ipix_right];
                    }
                }
            }
        }
    }

    // Convert the two spectra to one median.
    median(DIM_POL*ckd->dim_detector_spec,spectra.data(),strength);

    return 0;

} // }}}

int Swathcal::process_finalize( // {{{
)
{

    writelog(log_trace,"Evaluating ACT and ALT angles on L1B FOVs.");

    // These intermediate results are saved as members, purely as
    // optional detailed output.
    l1b_act_angles.resize(ckd->dim_fov); // Intermediate result.
    l1b_alt_angles.resize(ckd->dim_fov); // Intermediate result.
    // Running pointers for abscissa and actual angles sampled
    // on L1B FOVs.
    double *l1b_abscissa_cur = ckd->fov_act_angles.data();
    double *l1b_act_angles_cur = l1b_act_angles.data();
    double *l1b_alt_angles_cur = l1b_alt_angles.data();

    // Construct CKD.
    ckd->swath_swathvectors.resize(ckd->dim_fov*DIM_VEC,NC_FILL_DOUBLE);
    double *swathvector_cur = ckd->swath_swathvectors.data(); // Running pointer of the end result.
    ckd->swath_vectorplane_normals.resize(ckd->dim_vp*DIM_VEC,NC_FILL_DOUBLE);
    double *vectorplane_normal_cur = ckd->swath_vectorplane_normals.data(); // Running pointer of fitted plane.

    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {

        if (!ckd->vp_mask[ivp]) {

            // Number of data points, for all three of FOVcal ACT angles (abscissa),
            // swathcal ACT angles (ordinate) and swathcal ALT angles (ordinate).
            size_t &ndat = nbatch_vp[ivp];

            // Set mask and set larger weight for fits with higher magnitude.
            vector<double> local_s_y(ndat,1.0);
            vector<bool> local_mask(ndat); // Because you cannot offset vectors<bool>.
            size_t nfit_success = 0; // Meanwhile, count the number of successful fits.
            for (size_t idat=0 ; idat<ndat ; idat++) {
                local_mask[idat] = mask[ibatch_start[ivp]+idat];
                if (!local_mask[idat]) {
                    local_s_y[idat] = 1.0/pow(mag[ibatch_start[ivp]+idat],2.0);
                    nfit_success++;
                }
            }
            size_t order = 0; // Must be overwritten.
            for (size_t iattempt=0 ; iattempt<nattempt ; iattempt++) {
                if (nfit_success >= set->min_fit[iattempt]) {
                    order = set->order[iattempt];
                    break;
                }
            }
            check_error(order == 0,"Program Error: Impossible fit in viewport %zu should have been caught before.");
            writelog(log_debug,"Fitting ACT and ALT with %zu degrees of freedom for viewport %zu.",order,ivp);

            // Construct polynomial fitter.
            vector<double> knots = {-90.0*DEGREES,90*DEGREES}; // Knot values do not matter.
            Bspline b(this,order,2,knots.data());

            // Jacobian, common for ACT and ALT.
            vector<double> mat(b.nspline*ndat);
            handle(b.jaccalc(ndat,&act_fovcal[ibatch_start[ivp]],mat.data()));
            // Fit relationship for ACT angles.
            vector<double> res_act(b.nspline);
            check_error(linear_invert(b.nspline,ndat,mat.data(),OPT_DIAG,local_s_y.data(),&local_mask,&act_swathcal[ibatch_start[ivp]],res_act.data()) != 0,"Error: Failed fitting relationship between for ACT angles for viewport %zu. Possibly, some fit results were NaN.",ivp);
            // Evaluate the fit on the L1B FOVs.
            handle(b.jacapply(ckd->fov_nfov_vp[ivp],l1b_abscissa_cur,res_act.data(),l1b_act_angles_cur));

            // Repeat joke for ALT angles.
            vector<double> res_alt(b.nspline);
            check_error(linear_invert(b.nspline,ndat,mat.data(),OPT_DIAG,local_s_y.data(),&local_mask,&alt_swathcal[ibatch_start[ivp]],res_alt.data()) != 0,"Error: Failed fitting relationship between for ALT angles for viewport %zu. Possibly, some fit results were NaN.",ivp);
            // Evaluate the fit on the L1B FOVs.
            handle(b.jacapply(ckd->fov_nfov_vp[ivp],l1b_abscissa_cur,res_alt.data(),l1b_alt_angles_cur));

            // Turn these angles into swath vectors.
            for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {
                Vector res_irf; // Instrument reference frame.
                handle(theo->calculate_swath_vector(l1b_act_angles_cur[ifov],l1b_alt_angles_cur[ifov],res_irf));
                Vector res_srf = q_inst.quaternion_rotate(res_irf); // Satellite reference frame.
                res_srf.get(&swathvector_cur[ifov*DIM_VEC]); // Copy result to CKD.
            }

            // Fit a plane through all the swath vectors.
            // It is a plane that goes through (0,0,0).
            // In the ideal case, the equation is ax + cz = 0 and
            // y does not matter. In reality, y may matter a bit.
            // Let us use the model x + by + cz = 0, fixing the x
            // contribution to one so that the solution becomes unique.
            // This reasoning depends on the fact that the swath vectors
            // are already in satellite coordinates and that x is the
            // flight direction, z is nadir and y is to the right.
            // With just the instrument reference frame, this had not
            // been safe.

            // There are two state parameters and nfov_vp measurements.
            vector<double> jac(2*ckd->fov_nfov_vp[ivp]);
            vector<double> meas(ckd->fov_nfov_vp[ivp]);
            vector<double> res_bc(2);
            // The linear model will be: by + cz = -x with state vector
            // [b,c] and known parameters x,y,z.
            for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {
                double *xyz = &swathvector_cur[ifov*DIM_VEC]; // Swath vector as double pointer.
                jac[ifov] = xyz[1]; // d/db = y
                jac[ckd->fov_nfov_vp[ivp]+ifov] = xyz[2]; // d/dc = z
                meas[ifov] = -xyz[0]; // Measurement is -x.
            }
            check_error(linear_invert(2,ckd->fov_nfov_vp[ivp],jac.data(),OPT_NONE,NULL,NULL,meas.data(),res_bc.data()) != 0,"Error: Plane fit through swath vectors failed.")

            // Normalize, using the Vector class.
            Vector normal_cur(1.0,res_bc[0],res_bc[1]);
            normal_cur.normalize();

            // Write down results.
            normal_cur.get(vectorplane_normal_cur);

        }

        // Progress pointers.
        l1b_abscissa_cur += ckd->fov_nfov_vp[ivp];
        l1b_act_angles_cur += ckd->fov_nfov_vp[ivp];
        l1b_alt_angles_cur += ckd->fov_nfov_vp[ivp];
        swathvector_cur += ckd->fov_nfov_vp[ivp]*DIM_VEC;
        vectorplane_normal_cur += DIM_VEC;

    }

    return 0;

} // }}}

int Swathcal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{
    // Convert data type.
    vector<uint32_t> ifit_start_uint(ckd->dim_vp);
    vector<uint32_t> nfit_vp_uint(ckd->dim_vp);
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        ifit_start_uint[ivp] = ibatch_start[ivp];
        nfit_vp_uint[ivp] = nbatch_vp[ivp];
    }

    NcDim dimid_nfit;
    netcdf_check(nc,dimid_nfit = grp.addDim("fits",nbatch));
    netcdf_check(nc,grp.addVar("ifit_start",ncUint,ckd->dimid_vp).putVar(ifit_start_uint.data()));
    netcdf_check(nc,grp.addVar("nfit_vp",ncUint,ckd->dimid_vp).putVar(nfit_vp_uint.data()));
    netcdf_check(nc,grp.addVar("raw_act_fovcal",ncDouble,dimid_nfit).putVar(act_fovcal.data()));
    netcdf_check(nc,grp.addVar("raw_act_swathcal",ncDouble,dimid_nfit).putVar(act_swathcal.data()));
    netcdf_check(nc,grp.addVar("raw_alt_swathcal",ncDouble,dimid_nfit).putVar(alt_swathcal.data()));
    netcdf_check(nc,grp.addVar("raw_mag",ncDouble,dimid_nfit).putVar(mag.data()));
    netcdf_check(nc,grp.addVar("l1b_act_fovcal",ncDouble,ckd->dimid_fov).putVar(ckd->fov_act_angles.data()));
    netcdf_check(nc,grp.addVar("l1b_act_swathcal",ncDouble,ckd->dimid_fov).putVar(l1b_act_angles.data()));
    netcdf_check(nc,grp.addVar("l1b_alt_swathcal",ncDouble,ckd->dimid_fov).putVar(l1b_alt_angles.data()));

    return 0;

} // }}}

