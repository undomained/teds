// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "lininv.h"
#include "bspline.h"
#include "splitgaussfit.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "processor.h"
#include "wavecal.h"

// Settings functions.
Settings_wavecal::Settings_wavecal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "wave";
    opt.nonlin = false;
    opt.stray = false;
} // }}}
Settings_wavecal::~Settings_wavecal() {}
int Settings_wavecal::init_step( // {{{
    stringstream &stream,
    string &key,
    string &value,
    bool &recognized
)
{

    // Recognize specific settings.
    recognize_setting(peak_insignificant_island); // Setting for peak detection: Size of peaks or non-peaks to be discarded.
    recognize_setting(peak_maxreso); // Highest denominator in setting threshold to detect peaks.
    recognize_setting(peak_backgroundorder); // Degrees of freedom of background polynomial for peak detection.
    recognize_setting(peak_gaussfit); // Flag for performing gauss fit.
    recognize_setting(peak_domain_range); // Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    recognize_setting(peak_maxiter); // Number of iterations for fitting Gauss peak shape.
    recognize_setting(peak_magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
    recognize_setting(peak_postol); // Tolerance of peak positions (in indices).
    recognize_setting(peak_widthtol); // Tolerance of peak FWHM (in indices).
    recognize_setting(peak_backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
    recognize_setting(peak_steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
    recognize_setting(peak_initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
    recognize_setting(peak_reducerfactor_fail); // Multiplication of step reducer for a failed step.
    recognize_setting(peak_reducerfactor_success); // Multiplication of step reducer for a successful step.
    recognize_setting(peak_reducer_limit); // Maximum step control reducer for convergence.
    recognize_setting(refspec_insignificant_island); // In reference spectrum: Setting for peak detection: Size of peaks or non-peaks to be discarded.
    recognize_setting(refspec_maxreso); // In reference spectrum: Highest denominator in setting threshold to detect peaks.
    recognize_setting(refspec_backgroundorder); // In reference spectrum: Degrees of freedom of background polynomial for peak detection.
    recognize_setting(refspec_gaussfit); // In reference spectrum: Flag for performing gauss fit.
    recognize_setting(refspec_domain_range); // In reference spectrum: Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    recognize_setting(refspec_maxiter); // In reference spectrum: Number of iterations for fitting Gauss peak shape.
    recognize_setting(refspec_magtol); // In reference spectrum: Tolerance of peak magnitude (in units of a-priori magnitude).
    recognize_setting(refspec_postol); // In reference spectrum: Tolerance of peak positions (in indices).
    recognize_setting(refspec_widthtol); // In reference spectrum: Tolerance of peak FWHM (in indices).
    recognize_setting(refspec_backgroundtol); // In reference spectrum: Tolerance of background B-spline coefficients (in median signals).
    recognize_setting(refspec_steptolerance); // In reference spectrum: Inversion step acceptance criterion for fitting Gauss peak shape.
    recognize_setting(refspec_initial_step_reducer); // In reference spectrum: Starting value of the Levenberg-Marquardt step control parameter.
    recognize_setting(refspec_reducerfactor_fail); // In reference spectrum: Multiplication of step reducer for a failed step.
    recognize_setting(refspec_reducerfactor_success); // In reference spectrum: Multiplication of step reducer for a successful step.
    recognize_setting(refspec_reducer_limit); // In reference spectrum: Maximum step control reducer for convergence.
    recognize_setting(order); // Degrees of freedom in fit of wavelength against pixel index.
    recognize_setting(npeak); // Number of peak wavelengths.

    return 0;

} // }}}

// Constructor for the Wavecal CKD structure.
Wavecal::Wavecal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("wave");
    set = make_unique<Settings_wavecal>(this);
    Processor::set = set.get();
} // }}}
Wavecal::~Wavecal() {}

// 1. Put all the fitted peaks in a map.
// 2. Calibrate wavelengths in extracted spectrum perspective.
// 3. Define target wavelength grid.
int Wavecal::process_init( // {{{
)
{

    // Check legal settings.
    check_error(set->peak_maxreso < 2,"Error: Maximum peak-detection threshold resolution must be at least two. The higher this value (only powers of two make sense), the more patience the program has to find peaks. Setting 'peak_maxreso'");
    check_error(set->peak_backgroundorder < 1,"Error: B-spline degrees of freedom for the background for the Gauss fit must be at least one, for a constant background. Setting 'peak_backgroundorder'");
    // More inversion settings. Only relevant for a gauss fit.
    if (set->peak_gaussfit) {
        check_error(set->peak_domain_range == NC_FILL_DOUBLE,"Error: Missing domain range to determine fitting domain for the Guass, necessary because gaussfit is turned on. Setting 'peak_domain_range'");
        check_error(set->peak_maxiter < 1,"Error: Missing or zero maximum iterations in non-linear inversion. Setting 'peak_maxiter'");
        check_error(set->peak_magtol == NC_FILL_DOUBLE,"Error: Missing Gauss magnitude tolerance for convergence. Setting 'peak_magtol'");
        check_error(set->peak_postol == NC_FILL_DOUBLE,"Error: Missing Gauss position tolerance for convergence. Setting 'peak_postol'");
        check_error(set->peak_widthtol == NC_FILL_DOUBLE,"Error: Missing Gauss FWHM tolerance for convergence. Setting 'peak_widthtol'");
        check_error(set->peak_backgroundtol == NC_FILL_DOUBLE,"Error: Missing background tolerance in Gauss fit. Setting 'peak_backgroundtol'");
        check_error(set->peak_steptolerance < 1.0,"Error: Inversion step tolerance below one may lead to an endless loop. Setting 'peak_steptolerance'");
        check_error(set->peak_initial_step_reducer == NC_FILL_DOUBLE,"Error: Missing initial Levenber-Marquardt step reducer. Setting 'peak_initial_step_reducer'");
        if (set->peak_initial_step_reducer != 0.0) {
            check_error(set->peak_reducerfactor_success == NC_FILL_DOUBLE,"Error: Missing step reducer factor for successful steps. Settings 'peak_reducerfactor_success'");
            check_error(set->peak_reducerfactor_success > 1.0,"Error: Step reduction step for successful steps should be not higher than one. Settings 'peak_reducerfactor_success'");
        }
        if (set->peak_steptolerance == NC_FILL_DOUBLE) {
            // Step failure is impossible.
            if (set->peak_reducerfactor_fail != NC_FILL_DOUBLE) writelog(log_warning,"Warning: Step reducer factor for failed steps ignored, because no step tolerance is given, so steps cannot be rejected.");
        } else {
            // Check legal step reduction.
            check_error(set->peak_reducerfactor_fail == NC_FILL_DOUBLE,"Error: Missing reduction factor for failed steps. Stting 'peak_reducerfactor_fail'");
            check_error(set->peak_reducerfactor_fail <= 1.0,"Error: Step reduction factor for failed steps should be higher than one. Stting 'peak_reducerfactor_fail'");
        }
    }
    check_error(set->refspec_maxreso < 2,"Error: Maximum peak-detection threshold resolution must be at least two. The higher this value (only powers of two make sense), the more patience the program has to find peaks. Setting 'refspec_maxreso'");
    check_error(set->refspec_backgroundorder < 1,"Error: B-spline degrees of freedom for the background for the Gauss fit must be at least one, for a constant background. Setting 'refspec_backgroundorder'");
    // More inversion settings. Only relevant for a gauss fit.
    if (set->refspec_gaussfit) {
        check_error(set->refspec_domain_range == NC_FILL_DOUBLE,"Error: Missing domain range to determine fitting domain for the Guass, necessary because gaussfit is turned on. Setting 'refspec_domain_range'");
        check_error(set->refspec_maxiter < 1,"Error: Missing or zero maximum iterations in non-linear inversion. Setting 'refspec_maxiter'");
        check_error(set->refspec_magtol == NC_FILL_DOUBLE,"Error: Missing Gauss magnitude tolerance for convergence. Setting 'refspec_magtol'");
        check_error(set->refspec_postol == NC_FILL_DOUBLE,"Error: Missing Gauss position tolerance for convergence. Setting 'refspec_postol'");
        check_error(set->refspec_widthtol == NC_FILL_DOUBLE,"Error: Missing Gauss FWHM tolerance for convergence. Setting 'refspec_widthtol'");
        check_error(set->refspec_backgroundtol == NC_FILL_DOUBLE,"Error: Missing background tolerance in Gauss fit. Setting 'refspec_backgroundtol'");
        check_error(set->refspec_steptolerance < 1.0,"Error: Inversion step tolerance below one may lead to an endless loop. Setting 'refspec_steptolerance'");
        check_error(set->refspec_initial_step_reducer == NC_FILL_DOUBLE,"Error: Missing initial Levenber-Marquardt step reducer. Setting 'refspec_initial_step_reducer'");
        if (set->refspec_initial_step_reducer != 0.0) {
            check_error(set->refspec_reducerfactor_success == NC_FILL_DOUBLE,"Error: Missing step reducer factor for successful steps. Settings 'refspec_reducerfactor_success'");
            check_error(set->refspec_reducerfactor_success > 1.0,"Error: Step reduction step for successful steps should be not higher than one. Settings 'refspec_reducerfactor_success'");
        }
        if (set->refspec_steptolerance == NC_FILL_DOUBLE) {
            // Step failure is impossible.
            if (set->refspec_reducerfactor_fail != NC_FILL_DOUBLE) writelog(log_warning,"Warning: Step reducer factor for failed steps ignored, because no step tolerance is given, so steps cannot be rejected.");
        } else {
            // Check legal step reduction.
            check_error(set->refspec_reducerfactor_fail == NC_FILL_DOUBLE,"Error: Missing reduction factor for failed steps. Stting 'refspec_reducerfactor_fail'");
            check_error(set->refspec_reducerfactor_fail <= 1.0,"Error: Step reduction factor for failed steps should be higher than one. Stting 'refspec_reducerfactor_fail'");
        }
    }
    check_error(set->order < 1,"Error: B-spline degrees of freedom for curve over the detector must be at least one, for a constant spatial index per spectrum. Setting 'order'");
    check_error(set->npeak == 0,"Error: Missing (positive) number of peaks to recognize from image and reference spectrum.");

    // Shape the CKD.
    ckd->wave_spectra.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec,NC_FILL_DOUBLE);
    ckd->wave_target.resize(ckd->dim_fov*ckd->dim_detector_spec,NC_FILL_DOUBLE);

    // Peak fitter for spectra.
    spl = make_unique<Splitgaussfit>(this,set->npeak,set->peak_backgroundorder,ckd->dim_detector_spec);
    spl->setGaussfit(set->peak_gaussfit); // Flag for performing gauss fit.
    spl->setMaxiter(set->peak_maxiter); // Number of iterations for fitting Gauss peak shape.
    spl->setMagtol(set->peak_magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
    spl->setPostol(set->peak_postol); // Tolerance of peak positions (in indices).
    spl->setWidthtol(set->peak_widthtol); // Tolerance of peak FWHM (in indices).
    spl->setBackgroundtol(set->peak_backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
    spl->setSteptolerance(set->peak_steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
    spl->setInitialStepReducer(set->peak_initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
    spl->setReducerfactorFail(set->peak_reducerfactor_fail); // Multiplication of step reducer for a failed step.
    spl->setReducerfactorSuccess(set->peak_reducerfactor_success); // Multiplication of step reducer for a successful step.
    spl->setReducerLimit(set->peak_reducer_limit); // Maximum step control reducer for convergence.

    // Function fit of wavelenght as function of spectral index.
    vector<double> knots = {0.0,(double)ckd->dim_detector_spec-1};
    b = make_unique<Bspline>(this,set->order,2,knots.data());
    evalmat.resize(ckd->dim_detector_spec*b->nspline);
    vector<double> evalispec(ckd->dim_detector_spec);
    for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) evalispec[ispec] = (double)ispec;
    handle(b->jaccalc(ckd->dim_detector_spec,evalispec.data(),evalmat.data()));

    // Set viewport batches.
    ivp_batch.resize(ckd->dim_vp); // Maximum size.
    handle(batch_viewport(ivp_batch.data()));
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        // Wavelength calibration only works with one image per viewport.
        check_error(batches[ibatch].nl1a != 1,"Error: Wavelength calibration needs exactly one image per viewport. Viewport %zu has %zu.\n",ivp_batch[ibatch],batches[ibatch].nl1a);
    }

    return 0;

} // }}}
int Wavecal::process_batch( // {{{
    size_t ibatch // This is the viewport index.
)
{

    // The viewport index is linked to the batch index.
    size_t &ivp = ivp_batch[ibatch];

    // 1. Put all the fitted peaks in a map. {{{
    writelog(log_trace,"Fit peaks and put them into a map.");
    size_t npeak = set->npeak;
    vector<double> map_wavelength(ckd->fov_nfov_vp[ivp]*DIM_POL*npeak);
    vector<double> map_ispec(ckd->fov_nfov_vp[ivp]*DIM_POL*npeak);

    // Running pointers.
    double *ispec_cur = map_ispec.data(); // These are the positions that are fit.
    double *wave_cur = map_wavelength.data();

    L1A *l1a = l1a_instances[0];
    // Recognize peaks from reference spectrum.
    // In principle, any L1A file can have a differently-shaped reference spectrum.
    // So, we do not recycle any of the properties of this fitter.
    Splitgaussfit refspl(this,npeak,set->refspec_backgroundorder,l1a->dim_refspec);
    refspl.setX(l1a->refspec_wavelength.data());
    refspl.setGaussfit(set->refspec_gaussfit); // Flag for performing gauss fit.
    refspl.setMaxiter(set->refspec_maxiter); // Number of iterations for fitting Gauss peak shape.
    refspl.setMagtol(set->refspec_magtol); // Tolerance of peak magnitude (in units of a-priori magnitude).
    refspl.setPostol(set->refspec_postol); // Tolerance of peak positions (in indices).
    refspl.setWidthtol(set->refspec_widthtol); // Tolerance of peak FWHM (in indices).
    refspl.setBackgroundtol(set->refspec_backgroundtol); // Tolerance of background B-spline coefficients (in median signals).
    refspl.setSteptolerance(set->refspec_steptolerance); // Inversion step acceptance criterion for fitting Gauss peak shape.
    refspl.setInitialStepReducer(set->refspec_initial_step_reducer); // Starting value of the Levenberg-Marquardt step control parameter.
    refspl.setReducerfactorFail(set->refspec_reducerfactor_fail); // Multiplication of step reducer for a failed step.
    refspl.setReducerfactorSuccess(set->refspec_reducerfactor_success); // Multiplication of step reducer for a successful step.
    refspl.setReducerLimit(set->refspec_reducer_limit); // Maximum step control reducer for convergence.
    refspl.setMeas(l1a->refspec_radiance.data());
    // Reference spectrum has no noise or pixel mask.
    vector<double> res(3*npeak);
    check_error(refspl.solve(set->refspec_insignificant_island,set->refspec_domain_range,set->refspec_maxreso,res.data()) != 0,"Error: Unable to find the right number of peaks in reference spectrum.");
    vector<double> peaks(npeak);
    for (size_t ipeak=0 ; ipeak<npeak ; ipeak++) peaks[ipeak] = res[3*ipeak+1];
    for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {
        Spectra specs;
        handle(l1a->extract(ifov,specs));

        for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
            // Extract mask, because vector<bool> does not convert properly to
            // a bool pointer.
            vector<bool> mask_extracted(ckd->dim_detector_spec);
            for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) mask_extracted[ispec] = specs.mask[ipol*ckd->dim_detector_spec+ispec];
            fill_holes(ckd->dim_detector_spec,mask_extracted,&specs.signal[ipol*ckd->dim_detector_spec]);
            spl->setMeas(&specs.signal[ipol*ckd->dim_detector_spec]);
            spl->setNoise(&specs.noise[ipol*ckd->dim_detector_spec]);
            spl->setMask(&mask_extracted);
            vector<double> res(3*npeak);
            check_error(spl->solve(set->peak_insignificant_island,set->peak_domain_range,set->peak_maxreso,res.data()) != 0,"Error: Unable to find the right number of peaks in calibration measurement.");
            for (size_t ipeak=0 ; ipeak<npeak ; ipeak++) {
                ispec_cur[ipeak] = res[3*ipeak+1];
                wave_cur[ipeak] = peaks[ipeak];
            }
            // Move the pointers.
            ispec_cur += npeak;
            wave_cur += npeak;
        }
    }

    // }}}
    // 2. Calibrate wavelengths in extracted spectrum perspective. {{{
    writelog(log_trace,"Fit wavelengths per L1B field of view.");
    // For this, we fit a polynomial in ispec through all the fitted (ispec,wave)
    // points.
    // Reset moving pointers.
    ispec_cur = map_ispec.data();
    wave_cur = map_wavelength.data();
    // CKD pointer.
    double *wave_spectra_cur = &ckd->wave_spectra[ckd->dim_detector_spec];
    for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {
        for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
            vector<double> mat(npeak*b->nspline);
            handle(b->jaccalc(npeak,ispec_cur,mat.data()));
            vector<double> fit(b->nspline);
            check_error(linear_invert(b->nspline,npeak,mat.data(),OPT_NONE,NULL,NULL,wave_cur,fit.data()) != 0,"Error: Polynomial fit through fitted peaks failed.");
            // Manual matrix multiplication. The evaluation matrix is always the same,
            // so it is saved from viewport to viewport.
            for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                wave_spectra_cur[ispec] = 0.0;
                for (size_t ispline=0 ; ispline<b->nspline ; ispline++) {
                    wave_spectra_cur[ispec] += fit[ispline]*evalmat[ispline*ckd->dim_detector_spec+ispec];
                }
            }
            // Move the pointers.
            wave_spectra_cur += ckd->dim_detector_spec;
            ispec_cur += npeak;
            wave_cur += npeak;
        }
    }
    // }}}

    // 3. Define target wavelengths for co-sampling S+ and S-.
    // Reset spectra pointer.
    wave_spectra_cur = &ckd->wave_spectra[ckd->dim_detector_spec];
    // Set target-wavelength pointer.
    double *wave_target_cur = &ckd->wave_target[ckd->dim_detector_spec];
    for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {
        for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
            wave_target_cur[ispec] = 0.5*(wave_spectra_cur[ispec] + wave_spectra_cur[ckd->dim_detector_spec+ispec]);
        }
        // Move the pointers.
        wave_spectra_cur += DIM_POL*ckd->dim_detector_spec;
        wave_target_cur += ckd->dim_detector_spec;
    }

    return 0;

} // }}}

