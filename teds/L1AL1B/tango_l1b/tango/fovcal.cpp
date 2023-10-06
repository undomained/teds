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
#include "splitgaussfit.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "fovcal.h"

// Settings functions.
Settings_fovcal::Settings_fovcal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "fov";
    output_level_max = 1;
    opt.nonlin_apply = false;
} // }}}
Settings_fovcal::~Settings_fovcal() {}
int Settings_fovcal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(order); // Polynomial DOFs to fit the curves on the detector.
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
    recognize_setting(ispec_start); // First spectral pixel index to use for fitting polynomial.
    recognize_setting(ispec_end); // First spectral pixel index not to use for fitting polynomial.
    recognize_setting(outlier_cutoff); // As long as the difference between polynomial and data is this many times sigma, we remove the worst point.
    recognize_setting(hitthreshold); // Threshold for spectrum median intensity to be considered detected (fraction of highest).
    recognize_setting_vector(act_sampling); // ACT angle sampling per viewport.
    recognize_setting_vector(force_actangle_min); // Force minimum ACT angle as L1B FOV. Warning if not hit. This is one element per viewport.
    recognize_setting_vector(force_actangle_max); // Force maximum ACT angle as L1B FOV. Warning if not hit. This is one element per viewport.

    return 0;

} // }}}

// Constructor for the Fovcal CKD structure.
Fovcal::Fovcal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("fov");
    set = make_unique<Settings_fovcal>(this);
    Processor::set = set.get();
} // }}}
Fovcal::~Fovcal() {}

// Field-of-view calibration software.
// This is the planned protocol.
// 1. Extract the relevant metadata from the files.
// 2. Determine additional CKD dimensions.
// 3. Prepare the peak fitter.
// 4. Process per view port.
int Fovcal::process_init( // {{{
)
{

    // Check legal sttings. {{{
    check_error(set->order < 1,"Error: B-spline degrees of freedom for curve over the detector must be at least one, for a constant spatial index per spectrum. Setting 'order'");
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
    // Default domain end should be the end of the detector.
    if (set->ispec_end == NC_FILL_UINT64) set->ispec_end = ckd->dim_detector_spec; // Default is the entire detector.
    // Check fitting domain.
    check_error(set->ispec_end > ckd->dim_detector_spec,"Error: Spectral pixels outside the detector are selected to be taken into account. Setting 'ispec_end'");
    check_error(set->ispec_end <= set->ispec_start,"Error: No pixels selected for the fit. Make sure that the selected domain ends later than it begins. Settings 'ispec_start' and 'ispec_end'");
    // Threshold for hitting.
    check_error(set->hitthreshold >= 1.0,"Error: Threshold for a signal to be a hit should be less than one. Setting 'hitthreshold'");
    check_error(set->act_sampling.size() != ckd->dim_vp,"Error: Number of ACT angle samplings should be equal to the number of viewports. Setting 'act_sampling'");
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        check_error(set->act_sampling[ivp] <= 0.0,"Error: Zero or negative ACT angle sampling for viewport %zu. Setting 'act_sampling'",ivp);
        set->act_sampling[ivp] *= DEGREES; // Convert units.
    }
    // }}}

    // Peak fitter.
    spl = make_unique<Splitgaussfit>(this,2,set->backgroundorder,ckd->dim_detector_spat); // Two peaks.
    spl->setGaussfit(set->gaussfit); // Flag for performing gauss fit.
    spl->setMaxiter(set->maxiter); // Number of iterations for fitting Gauss peak shape.
    spl->setMagtol(set->magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
    spl->setPostol(set->postol); // Tolerance of peak positions (in indices).
    spl->setWidthtol(set->widthtol); // Tolerance of peak FWHM (in indices).
    spl->setBackgroundtol(set->backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
    spl->setSteptolerance(set->steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
    spl->setInitialStepReducer(set->initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
    spl->setReducerfactorFail(set->reducerfactor_fail); // Multiplication of step reducer for a failed step.
    spl->setReducerfactorSuccess(set->reducerfactor_success); // Multiplication of step reducer for a successful step.
    spl->setReducerLimit(set->reducer_limit); // Maximum step control reducer for convergence.

    vector<double> knots = {0.0,(double)ckd->dim_detector_spec-1};
    b = make_unique<Bspline>(this,set->order,2,knots.data());
    mat_full.resize(b->nspline*ckd->dim_detector_spec);
    vector<double> ispec_full(ckd->dim_detector_spec);
    for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) ispec_full[ispec] = (double)ispec;
    handle(b->jaccalc(ckd->dim_detector_spec,ispec_full.data(),mat_full.data()));

    // Acquire maximum number of FOV angles.
    // This assumes that all measurements are visible.
    // Meanwhile, check that all L1A files are from one viewport only.
    size_t cnt = 0;
    uint8_t l1a_viewport_cur = 1;
    size_t nfov_max = 0;
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        double mn_act = NC_FILL_DOUBLE;
        double mx_act = NC_FILL_DOUBLE;
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            L1A *l1a = l1a_instances_all[il1a].get();
            if (l1a->viewport != l1a_viewport_cur) continue;
            cnt++;
            double &act_cur = l1a->act_angle;
            if (mn_act == NC_FILL_DOUBLE || act_cur < mn_act) mn_act = act_cur;
            if (mx_act == NC_FILL_DOUBLE || act_cur > mx_act) mx_act = act_cur;
        }
        nfov_max += static_cast<size_t>((mx_act-mn_act)/set->act_sampling[ivp]) + 1;
        l1a_viewport_cur *= 2;
    }
    check_error(cnt != nl1a_total,"Error: Not all L1A files are from exactly one viewport.");

    // Shape fixed-sized CKD.
    ckd->fov_nfov_vp.resize(ckd->dim_vp,0);
    ckd->dim_fov = 0; // This is a counter.
    // Organizational arrays. These have fixed contents.
    ckd->fov_dims_spec.resize(nfov_max,ckd->dim_detector_spec); // Fixed known contents for unbinned CKD.
    ckd->fov_iel_start.resize(nfov_max);
    for (size_t ifov=0 ; ifov<nfov_max ; ifov++) ckd->fov_iel_start[ifov] = ifov*DIM_POL*ckd->dim_detector_spec;
    // Allocate the maximum size for the CKD.
    ckd->fov_act_angles.resize(nfov_max);
    ckd->fov_ispat.resize(nfov_max*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_ipix1.resize(nfov_max*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_ipix2.resize(nfov_max*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_weight1.resize(nfov_max*DIM_POL*ckd->dim_detector_spec);

    // Set viewport batches.
    ivp_batch.resize(ckd->dim_vp); // Maximum size.
    handle(batch_viewport(ivp_batch.data()));
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        check_error(batches[ibatch].nl1a == 0,"Error: Viewport %zu has no input L1A images.",ivp_batch[ibatch]);
    }

    if (set->output_level >= 1) {
        detailed_output_viewports.resize(ckd->dim_vp);
        for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
            // The viewport index is linked to the batch index.
            size_t &ivp = ivp_batch[ibatch];
            // The output instances will be indexed with viewport indices.
            // This means that some empty objects may exist.
            // It is not default to skip viewports, and the index ivp is
            // more logical than ibatch.
            DetailedOutputFov &det_out = detailed_output_viewports[ivp];
            size_t &nl1a_cur = batches[ibatch].nl1a;
            det_out.nl1a = nl1a_cur;
            det_out.act_angles.resize(nl1a_cur); // ACT angles of the files.
            det_out.ispat_raw.resize(nl1a_cur*DIM_POL*ckd->dim_detector_spec); // Raw fit results. Dimension (L1A,pol,detector_spec).
            det_out.ispat_smooth.resize(nl1a_cur*DIM_POL*ckd->dim_detector_spec); // Smoothed fit results, still per L1A file. Dimension (L1A,pol,detector_spec).
            det_out.spectrummedians.resize(nl1a_cur); // Median spectra for visibility assessment. Dimension (L1A).
            det_out.visible.resize(nl1a_cur); // Visibility flag (L1A).
        }
    }

    return 0;

} // }}}

int Fovcal::process_batch(size_t ibatch, const Calibration_options& opt)
{

    // Set up administration. {{{

    // The viewport index is linked to the batch index.
    size_t &ivp = ivp_batch[ibatch];

    // Evaluated fit results on measured angles.
    vector<double> meas_ispat(nl1a*DIM_POL*ckd->dim_detector_spec);

    // Running pointers.
    double *meas_ispatm_cur = meas_ispat.data();
    double *meas_ispatp_cur = meas_ispat.data() + ckd->dim_detector_spec;

    // Detailed output pointer.
    // This is the one that filled within the loop.
    double *det_ispat_raw = NULL;
    if (set->output_level >= 1) det_ispat_raw = detailed_output_viewports[ivp].ispat_raw.data();

    // Sort L1A-images by their across-track angles. This is needed for
    // later interpolation on L1B ACT angles.
    vector<double> act_angles_unsorted(nl1a);
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) act_angles_unsorted[il1a] = l1a_instances[il1a]->act_angle;
    vector<double> act_angles_sorted(nl1a);
    vector<size_t> l1a_indexarray(nl1a);
    Array<double>::sort(nl1a,act_angles_unsorted.data(),act_angles_sorted.data(),l1a_indexarray.data());

    // Detailed output: Sorted ACT angles.

    // If the detector is missed, you see nothing or a very little bit. If
    // you see nothing, we expect a lot of non-fits and chaos. If we see a
    // a little bit, we may get a good fit with all the smoothing included,
    // but the signal is low. If you get chaos, the signal is also low.
    // So if we take the medians of the spectrum that we think that belongs
    // to the measurement and compare it to the maximum of it, we have
    // a measure. It would be a shame if there are missing FOV angles
    // between hitting ones. That is an error or a warning.
    vector<double> spectrummedians(nl1a);
    // These are the median of S+ and S- together.
    double *med_cur = spectrummedians.data(); // Running pointer.
    double spectrummedians_max = 0.0;

    // }}}

    percentagelog_open("Processing per calibration measurement");
    for (size_t il1a = 0 ; il1a<nl1a ; il1a++) {

        // Process per measurement. {{{
        percentagelog_progress(il1a,nl1a);

        // Set pointer to appropriate L1A instance.
        L1A *l1a = l1a_instances[l1a_indexarray[il1a]];

        // Retrieve peak positions. {{{

        // Fitted spatial pixel indices and their success flags.
        // By default everything fails. Only if the end of the calculation
        // is reached, the flag is set to success.
        vector<double> fit_spatm(ckd->dim_detector_spec);
        vector<bool> mask_spatm(ckd->dim_detector_spec,true);
        vector<double> fit_spatp(ckd->dim_detector_spec);
        vector<bool> mask_spatp(ckd->dim_detector_spec,true);
        // Both have their own mask. The mask may get polluted during the
        // outlier protection function.

        // Detailed output. Write fill values in the array of raw fits.
        // Non-processed spectral indices as well as failed fit will
        // stay fillvalue.
        if (set->output_level >= 1) {
            for (size_t iel=0 ; iel<DIM_POL*ckd->dim_detector_spec ; iel++) det_ispat_raw[iel] = NC_FILL_DOUBLE;
        }

        // Loop over lines with same spectral index.
        for (size_t ispec = set->ispec_start ; ispec < set->ispec_end ; ispec++) {

            // Check if any pixel lives. If no pixel lives, you die.
            bool living = false;
            for (size_t ipix=ispec ; ipix<ckd->npix ; ipix+=ckd->dim_detector_spec) {
                if (!l1a->pixelmask[ipix]) {
                    living = true;
                    break;
                }
            }
            if (!living) continue; // Die.

            // Extract line. {{{

            // Dead pixels are set to the average of their
            // living neighbours. There is at least one living neighbour
            // somewhere.
            vector<double> line(ckd->dim_detector_spat);
            vector<double> linenoise(ckd->dim_detector_spat);
            vector<bool> linemask(ckd->dim_detector_spat);
            for (size_t ispat=0 ; ispat<ckd->dim_detector_spat ; ispat++) {
                size_t ipix = ispat*ckd->dim_detector_spec+ispec;
                line[ispat] = l1a->image[ipix];
                linenoise[ispat] = l1a->noise[ipix];
                linemask[ispat] = l1a->pixelmask[ipix];
            }
            fill_holes(ckd->dim_detector_spat,linemask,line.data());
            // No need to fill noise holes. The fill_holes routine is needed to fill the holes
            // of the signal for the peak detection system. There, the noise and the mask cannot
            // be used, so holes have to be estimated. For all inversions, the mask is used, so
            // the noise at dead pixels does not matter.

            // }}}

            // Perform the retrieval of the peaks.
            // If any action fails, we just mark failure and continue
            // to the next pixel.
            spl->setMeas(line.data());
            spl->setNoise(linenoise.data(),false); // Flag is the S-flag. It is sigma, not S.
            spl->setMask(&linemask);
            vector<double> res(6);
            if (spl->solve(set->insignificant_island,set->domain_range,set->maxreso,res.data()) != 0) continue;

            // Filter that the result is inside the detector.
            if (res[1] < 0.0 || res[1] > ckd->dim_detector_spat-1.0 || res[4] < 0.0 || res[4] > ckd->dim_detector_spat-1.0) continue;

            // Write down result and mark success.
            fit_spatm[ispec] = res[1];
            mask_spatm[ispec] = false;
            fit_spatp[ispec] = res[4];
            mask_spatp[ispec] = false;
            if (set->output_level >= 1) {
                // Write raw result.
                det_ispat_raw[ispec] = res[1];
                det_ispat_raw[ckd->dim_detector_spec+ispec] = res[4];
            }
        }
        // }}}

        *med_cur = 0.0; // Use this value when the loop is broken.
        // Breakable do-once loop.
        for (bool ctd=true ; ctd ; ctd=false) {
            vector<double> spline_m(b->nspline);
            vector<double> spline_p(b->nspline);
            if (set->outlier_cutoff == NC_FILL_DOUBLE) {
                if (linear_invert(b->nspline,ckd->dim_detector_spec,mat_full.data(),OPT_NONE,NULL,&mask_spatm,fit_spatm.data(),spline_m.data())) break;
                if (linear_invert(b->nspline,ckd->dim_detector_spec,mat_full.data(),OPT_NONE,NULL,&mask_spatp,fit_spatp.data(),spline_p.data())) break;
            } else {
                if (linear_invert_outlierprotect(b->nspline,ckd->dim_detector_spec,mat_full.data(),OPT_NONE,NULL,&mask_spatm,fit_spatm.data(),set->outlier_cutoff,spline_m.data())) break;
                if (linear_invert_outlierprotect(b->nspline,ckd->dim_detector_spec,mat_full.data(),OPT_NONE,NULL,&mask_spatp,fit_spatp.data(),set->outlier_cutoff,spline_p.data())) break;
            }

            // Evaluate the fit.
            Matrix::matmul_rl_fold_slow(b->nspline,ckd->dim_detector_spec,mat_full.data(),spline_m.data(),meas_ispatm_cur);
            Matrix::matmul_rl_fold_slow(b->nspline,ckd->dim_detector_spec,mat_full.data(),spline_p.data(),meas_ispatp_cur);

            // Extract the spectra to get the median.
            vector<double> spectra(DIM_POL*ckd->dim_detector_spec);
            double *meas_ispat_cur = meas_ispatm_cur;
            for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                    double ispat = *meas_ispat_cur;
                    int ispat_left = static_cast<int>(ispat);
                    double weightleft = ispat_left+1-ispat;
                    size_t ipix_left = ispat_left*ckd->dim_detector_spec + ispec;
                    size_t ipix_right = (ispat_left+1)*ckd->dim_detector_spec + ispec;
                    if (ispat_left < 0 || ispat_left > (int)ckd->dim_detector_spat-2) {
                        spectra[ipol*ckd->dim_detector_spec+ispec] = 0.0;
                    } else {
                        if (l1a->pixelmask[ipix_left]) {
                            if (l1a->pixelmask[ipix_right]) {
                                spectra[ipol*ckd->dim_detector_spec+ispec] = 0.0;
                            } else {
                                spectra[ipol*ckd->dim_detector_spec+ispec] = l1a->image[ipix_right];
                            }
                        } else {
                            if (l1a->pixelmask[ipix_right]) {
                                spectra[ipol*ckd->dim_detector_spec+ispec] = l1a->image[ipix_left];
                            } else {
                                spectra[ipol*ckd->dim_detector_spec+ispec] = weightleft * l1a->image[ipix_left] + (1.0-weightleft) * l1a->image[ipix_right];
                            }
                        }
                    }
                    // Move the pointer.
                    meas_ispat_cur++;
                }
            }
            // Convert the two spectra to one median.
            median(DIM_POL*ckd->dim_detector_spec,spectra.data(),*med_cur);
            if (*med_cur > spectrummedians_max) spectrummedians_max = *med_cur;
        }

        // Progress the running pointers to the next measurement.
        meas_ispatm_cur += DIM_POL*ckd->dim_detector_spec;
        meas_ispatp_cur += DIM_POL*ckd->dim_detector_spec;
        med_cur++;
        if (set->output_level >= 1) det_ispat_raw += DIM_POL*ckd->dim_detector_spec;

        // }}}

    }
    percentagelog_close();

    // Detailed output: Write smoothed fits and medians.
    if (set->output_level >= 1) {
        DetailedOutputFov &det_out = detailed_output_viewports[ivp];
        // Copy (sorted) ACT angles.
        memcpy(det_out.act_angles.data(),act_angles_sorted.data(),nl1a*sizeof(double));
        // Copy spectrum medians.
        memcpy(det_out.spectrummedians.data(),spectrummedians.data(),nl1a*sizeof(double));
        // Copy smoothed spectra.
        memcpy(det_out.ispat_smooth.data(),meas_ispat.data(),nl1a*DIM_POL*ckd->dim_detector_spec*sizeof(double));
    }

    // The medians should be high in the middle and low on the sides.
    // The threshold of being high must be some percentage of the maximum
    // from the settings.
    // Find the islands of hit spectra.
    size_t il1a_start_largest_island = NC_FILL_UINT64;
    size_t largest_island = 0;
    size_t current_island = 0;
    size_t total = 0;
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
        if (spectrummedians[il1a] > set->hitthreshold*spectrummedians_max) {
            total++;
            current_island++;
            if (set->output_level >= 1) detailed_output_viewports[ivp].visible[il1a] = 1;
        } else {
            if (current_island > largest_island) {
                largest_island = current_island;
                il1a_start_largest_island = il1a-largest_island;
            }
            current_island = 0;
            if (set->output_level >= 1) detailed_output_viewports[ivp].visible[il1a] = 0;
        }
    }
    // Close off any last island if it extends to the end.
    if (current_island > largest_island) {
        largest_island = current_island;
        il1a_start_largest_island = nl1a-largest_island;
    }
    check_error(il1a_start_largest_island == NC_FILL_UINT64,"Program error: Variable il1a_start_largest_island not properly obtained during search algorithm.");
    if (total > largest_island) {
        writelog(log_warning,"Warning: Across-track angles that hit the detector seem not to be contiguous.");
    }
    check_error(total == 0,"Program error: No visible spectrum detected. This should never happen.");

    // Find out lowest and highest angles that hit the detetor.
    double act_angle_min = act_angles_sorted[il1a_start_largest_island];
    double act_angle_max = act_angles_sorted[il1a_start_largest_island+largest_island-1];

    // Overwrite boundaries by chosen force. Warn if the force is outside
    // the part that hits the detector. Note that you can still be between
    // the last hitting and the first missing, but that is still worth a
    // warning.
    if (set->force_actangle_min.size() > ivp) {
        double force = set->force_actangle_min[ivp];
        if (force >= -180.0 && force <= 180.0) {
            if (force*DEGREES < act_angle_min) {
                writelog(log_warning,"Warning: Forced lowest ACT angle : %.8f degrees does not seem to hit the detector. Minimum hitting ACT angle is %.8f degrees.",force,act_angle_min/DEGREES);
                act_angle_min = force*DEGREES;
            }
        }
    }
    if (set->force_actangle_max.size() > ivp) {
        double force = set->force_actangle_max[ivp];
        if (force >= -180.0 && force <= 180.0) {
            if (force*DEGREES > act_angle_max) {
                writelog(log_warning,"Warning: Forced highest ACT angle : %.8f degrees does not seem to hit the detector. Maximum hitting ACT angle is %.8f degrees.",force,act_angle_max/DEGREES);
                act_angle_max = force*DEGREES;
            }
        }
    }
    check_error(act_angle_max < act_angle_min,"Error: No ACT angle domain left for viewport %zu, probably through force. Derived range is from %.7f to %.7f.",ivp,act_angle_min/DEGREES,act_angle_max/DEGREES);

    // Define L1B swaths.
    uint32_t &nfov_cur = ckd->fov_nfov_vp[ivp];
    nfov_cur = static_cast<size_t>((act_angle_max-act_angle_min)/set->act_sampling[ivp]) + 1;
    double act_start = 0.5*(act_angle_min + act_angle_max - (nfov_cur-1)*set->act_sampling[ivp]);
    for (size_t ifov=0 ; ifov<nfov_cur ; ifov++) {
        ckd->fov_act_angles[ifov] = act_start + ifov*set->act_sampling[ivp];
    }
    ckd->dim_fov += nfov_cur;

    // Interpolate on L1B FOVs. {{{
    writelog(log_trace,"Interpolate on L1B FOVs.");

    // Running CKD pointers.
    double *fov_ispat_cur = &ckd->fov_ispat[ckd->dim_detector_spec];
    size_t *fov_ipix1_cur = &ckd->fov_ipix1[ckd->dim_detector_spec];
    size_t *fov_ipix2_cur = &ckd->fov_ipix2[ckd->dim_detector_spec];
    double *fov_weight1_cur = &ckd->fov_weight1[ckd->dim_detector_spec];

    // We will first get the indices and their weights and execute the
    // interpolation for each detector pixel.
    for (size_t ifov=0 ; ifov<nfov_cur ; ifov++) {
        if (largest_island == 1) {
            // This is no interpolation. Interpolation may lead to a
            // segmentation fault with zero times element -1.
            size_t iel = 0; // Folded iterator.
            for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                    double &ispat = meas_ispat[il1a_start_largest_island*DIM_POL*ckd->dim_detector_spec+iel];
                    fov_ispat_cur[iel] = ispat;
                    check_error(ispat < 0.0 || ispat > static_cast<double>(ckd->dim_detector_spec)-1.0,"Error: Spectrum positioning seems to be outside the detector, spatial pixel index %.7f",ispat);
                    // This is unbinned CKD, so the relationship between
                    // coordinates and pixel number is governed by the
                    // detector size.
                    fov_ipix1_cur[iel] = static_cast<size_t>(ispat)*ckd->dim_detector_spec + ispec;
                    fov_ipix2_cur[iel] = fov_ipix1_cur[iel] + ckd->dim_detector_spec;
                    fov_weight1_cur[iel] = floor(ispat) + 1.0 - ispat;
                    iel++; // Increase folded iterator.
                }
            }
        } else {
            size_t il1a_left = NC_FILL_UINT64;
            double ang = ckd->fov_act_angles[ifov];
            if (ang <= act_angles_sorted[il1a_start_largest_island]) il1a_left = il1a_start_largest_island;
            else if (ang >= act_angles_sorted[il1a_start_largest_island+largest_island-1]) il1a_left = il1a_start_largest_island+largest_island-2;
            else {
                for (size_t il1a_right = il1a_start_largest_island+1 ; il1a_right < il1a_start_largest_island+largest_island ; il1a_right++) {
                    if (act_angles_sorted[il1a_right] >= ang) {
                        il1a_left = il1a_right - 1;
                        break;
                    }
                }
            }
            check_error(il1a_left == NC_FILL_UINT64,"Program error: Variable il1a_left not properly initialized during search loop.");
            // Get the weight factor. This implies linear extrapolation if needed.
            // If that is needed, a warning should have been raised before.
            double weightleft = (act_angles_sorted[il1a_left+1]-ang) / (act_angles_sorted[il1a_left+1]-act_angles_sorted[il1a_left]);
            // Apply interpolation.
            size_t iel = 0; // Folded iterator.
            for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
                for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                    double ispat = weightleft * meas_ispat[il1a_left*DIM_POL*ckd->dim_detector_spec+iel] + (1.0-weightleft) * meas_ispat[(il1a_left+1)*DIM_POL*ckd->dim_detector_spec+iel];
                    check_error(ispat < 0.0 || ispat > static_cast<double>(ckd->dim_detector_spat)-1.0,"Error: Spectrum positioning seems to be outside the detector, spatial pixel index %.7f",ispat);
                    fov_ispat_cur[iel] = ispat;

                    // Convert iterpolated ispat to two pixels and a weight.
                    fov_ipix1_cur[iel] = static_cast<size_t>(ispat)*ckd->dim_detector_spec + ispec;
                    fov_ipix2_cur[iel] = fov_ipix1_cur[iel] + ckd->dim_detector_spec;
                    fov_weight1_cur[iel] = floor(ispat) + 1.0 - ispat;
                    iel++; // Increase folded iterator.
                }
            }
        }
        // Progress running pointers.
        fov_ispat_cur += DIM_POL*ckd->dim_detector_spec;
        fov_ipix1_cur += DIM_POL*ckd->dim_detector_spec;
        fov_ipix2_cur += DIM_POL*ckd->dim_detector_spec;
        fov_weight1_cur += DIM_POL*ckd->dim_detector_spec;

    }
    // }}}

    return 0;

} // }}}

int Fovcal::process_finalize( // {{{
)
{
    // Clip the arrays to their right sizes. This is not necessary, but cleaner. If anyone ever
    // uses the size of these arrays, this clipping becomes necessary.
    ckd->fov_dims_spec.resize(ckd->dim_fov);
    ckd->fov_iel_start.resize(ckd->dim_fov);
    ckd->fov_act_angles.resize(ckd->dim_fov);
    ckd->fov_ispat.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_ipix1.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_ipix2.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec);
    ckd->fov_weight1.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec);

    return 0;

} // }}}

int Fovcal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{

    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) {
        if (ckd->vp_mask[ivp]) continue; // Skip any skipped viewport.
        DetailedOutputFov &det_out = detailed_output_viewports[ivp];
        string viewport_suffix = format("_%zu",ivp);
        string nm = format("nmeas%s",viewport_suffix.c_str());
        NcDim dimid_nmeas;
        netcdf_check(nc,dimid_nmeas = grp.addDim(nm.c_str(),det_out.nl1a));
        nm = format("act_angles%s",viewport_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dimid_nmeas).putVar(det_out.act_angles.data()));
        vector<NcDim> dims = {dimid_nmeas,ckd->dimid_pol,ckd->dimid_detector_spec};
        nm = format("ispat_raw%s",viewport_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dims).putVar(det_out.ispat_raw.data()));
        nm = format("ispat_smooth%s",viewport_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dims).putVar(det_out.ispat_smooth.data()));
        nm = format("spectrummedians%s",viewport_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncDouble,dimid_nmeas).putVar(det_out.spectrummedians.data()));
        nm = format("visible%s",viewport_suffix.c_str());
        netcdf_check(nc,grp.addVar(nm.c_str(),ncUbyte,dimid_nmeas).putVar(det_out.visible.data()));
    }

    return 0;

} // }}}

