// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "lininv.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "ckd.h"
#include "l1a.h"
#include "prnucal.h"

// Settings functions.
Settings_prnucal::Settings_prnucal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "prnu";
    output_level_max = 1;
} // }}}
Settings_prnucal::~Settings_prnucal() {}
int Settings_prnucal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    recognize_setting(order_spat); // Low-order polynomial DOFs for spatial dimension (one higher than order).
    recognize_setting(order_spec); // Low-order polynomial DOFs for spectral dimension (one higher than order).
    recognize_setting(outlier_cutoff); // Outlier-protection threshold for low-order polynomial fit.
    recognize_setting(iterative); // Flag for letting the pixel mask and the PRNU iteratively converge.
    recognize_setting(mask_prnu_min); // Minimum pixel response non-uniformity.
    recognize_setting(mask_prnu_max); // Maximum pixel response non-uniformity.

    return 0;

} // }}}

// Constructor for the Prnucal CKD structure.
Prnucal::Prnucal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("prnu");
    set = make_unique<Settings_prnucal>(this);
    Processor::set = set.get();
    own_skip = &ckd->prnu_skip;
} // }}}
Prnucal::~Prnucal() {}

int Prnucal::process_init( // {{{
)
{

    // Verify legal settings.
    check_error(set->order_spat < 1,"Error: B-spline degrees of freedom for smooth fit in spatial dimension must be at least one, for a constant fit. Setting 'order_spat'");
    check_error(set->order_spec < 1,"Error: B-spline degrees of freedom for smooth fit in spatial dimension must be at least one, for a constant fit. Setting 'order_spec'");

    // Verify existence of necessary GSE data.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        check_error(l1a_instances_all[il1a]->illumination_level == NC_FILL_DOUBLE,"Error: Missing illumination level as GSE-data from '%s'.",l1a_instances_all[il1a]->filename.c_str());
    }

    return 0;

} // }}}

// PRNU correction calibration software.
// This is the planned protocol.
int Prnucal::process_batch( // {{{
    size_t ibatch // Meaningless zero here.
)
{

    // Shape the output.
    ckd->prnu_prnu.resize(ckd->npix);

    // Save nl1a in a safe way for detailed output.
    if (set->output_level >= 1) det_nl1a = nl1a;

    // After correcting for non-linearity and the exposure time, the signal should
    // be linear in the number of photons (intensity times exposure time), but the
    // proportionality constant differs per pixel. That is the PRNU. We will calibrate
    // out the quick osciallation, but keep the broad ones.
    photons.resize(nl1a); // Save array for detailed output, otherwise, just use it.
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) photons[il1a] = l1a_instances[il1a]->illumination_level * l1a_instances[il1a]->exposure_time;

    // The sensitivity is the signal divided by the number of photons. The idea is
    // that the sensitivity is PRNU times large-scale sensitivity, where the large-
    // scale sensitivity is considered not to be pixel, but light related.
    sensitivity.resize(ckd->npix,0.0); // Save array for detailed output, otherwise, just use it.
    // The sensitivity is averaged over all the measurements.
    for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
        for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) sensitivity[ipix] += (l1a_instances[il1a]->image[ipix] / photons[il1a]) / nl1a;
    }

    // Because the detector grid is a rectangular grid, the polynomial fit can be
    // done by dimension as well. With many pixels, it seems better to do that.
    writelog(log_trace,"Constructing 1D Jacobians.");
    vector<double> knots_spat = {0.0,(double)ckd->dim_detector_spat-1.0};
    Bspline b_spat(this,set->order_spat,2,knots_spat.data());
    vector<double> jac_spat(b_spat.nspline*ckd->dim_detector_spat);
    vector<double> idx_spat(ckd->dim_detector_spat);
    for (size_t ispat=0 ; ispat<ckd->dim_detector_spat ; ispat++) idx_spat[ispat] = (double) ispat;
    handle(b_spat.jaccalc(ckd->dim_detector_spat,idx_spat.data(),jac_spat.data()));

    vector<double> knots_spec = {0.0,(double)ckd->dim_detector_spec-1.0};
    Bspline b_spec(this,set->order_spec,2,knots_spec.data());
    vector<double> jac_spec(b_spec.nspline*ckd->dim_detector_spec);
    vector<double> idx_spec(ckd->dim_detector_spec);
    for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) idx_spec[ispec] = (double) ispec;
    handle(b_spec.jaccalc(ckd->dim_detector_spec,idx_spec.data(),jac_spec.data()));

    // Construct aggregate Jacobian.
    writelog(log_trace,"Constructing aggregate Jacobians.");
    size_t nstate = b_spat.nspline * b_spec.nspline;
    vector<double> jac_2d(nstate*ckd->npix);
    double *jacptr = jac_2d.data();
    for (size_t ispline_spat=0 ; ispline_spat<b_spat.nspline ; ispline_spat++) {
        for (size_t ispline_spec=0 ; ispline_spec<b_spec.nspline ; ispline_spec++) {
            for (size_t ispat=0 ; ispat<ckd->dim_detector_spat ; ispat++) {
                for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                    *jacptr = jac_spat[ispline_spat*ckd->dim_detector_spat+ispat] * jac_spec[ispline_spec*ckd->dim_detector_spec+ispec];
                    jacptr++;
                }
            }
        }
    }
    // Optional iterative loop.
    sensitivity_fit.resize(ckd->npix); // Save array for detailed output, otherwise, just use it.
    for (bool ctd=true ; ctd ; ctd=set->iterative) {
        // Apply prior pixel mask. All masked pixels are excluded.
        vector<bool> excluded(ckd->npix);
        for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) excluded[ipix] = ckd->mask[ipix];
        writelog(log_trace,"Performing the fit.");
        vector<double> res(nstate);
        if (set->outlier_cutoff == NC_FILL_DOUBLE) {
            check_error(linear_invert(nstate,ckd->npix,jac_2d.data(),OPT_NONE,NULL,&excluded,sensitivity.data(),res.data()) != 0,"Error: Smooth fit of detector stimulus failed.");
        } else {
            check_error(linear_invert_outlierprotect(nstate,ckd->npix,jac_2d.data(),OPT_NONE,NULL,&excluded,sensitivity.data(),set->outlier_cutoff,res.data()) != 0,"Error: Smooth fit of detector stimulus failed.");
        }
        Matrix::matmul_rl_fold_slow(nstate,ckd->npix,jac_2d.data(),res.data(),sensitivity_fit.data());

        // The PRNU is defined as the sensitivity divided by the fitted (low-order)
        // sensitivity.
        writelog(log_trace,"Interpreting results.");
        bool converged = true;
        for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
            ckd->prnu_prnu[ipix] = sensitivity[ipix] / sensitivity_fit[ipix];
            // Apply mask.
            if (
                (set->mask_prnu_max != NC_FILL_DOUBLE && ckd->prnu_prnu[ipix] > set->mask_prnu_max) ||
                (set->mask_prnu_min != NC_FILL_DOUBLE && ckd->prnu_prnu[ipix] < set->mask_prnu_min)
            ) {
                if (!excluded[ipix]) converged = false;
                ckd->mask[ipix] = true;
            }
        }
        if (converged) break;
    }

    return 0;

} // }}}

int Prnucal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{
    NcDim dimid_nmeas;
    netcdf_check(nc,dimid_nmeas = grp.addDim("nmeas",det_nl1a));
    NcVar var_photons;
    netcdf_check(nc,var_photons = grp.addVar("photons",ncDouble,dimid_nmeas));
    netcdf_check(nc,var_photons.putAtt("meaning","Ilummination level multiplied by exposure time"));
    netcdf_check(nc,var_photons.putVar(photons.data()));
    vector<NcDim> dims = {ckd->dimid_detector_spat,ckd->dimid_detector_spec};
    netcdf_check(nc,grp.addVar("sensitivity",ncDouble,dims).putVar(sensitivity.data()));
    netcdf_check(nc,grp.addVar("sensitivity_fit",ncDouble,dims).putVar(sensitivity_fit.data()));

    return 0;

} // }}}

