#include "header.h"
#include "functions.h"
#include "fourier.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "settings_l1b.h"
#include "settings_isrf.h"
#include "settings_geo.h"
#include "planet.h"
#include "utc.h"
#include "ckd.h"
#include "l1x_inputfile.h"
#include "frame.h"
#include "vector.h"
#include "linear_spline.h"
#include "cubic_spline.h"
#include <numbers>

#include <numeric>
#include <cmath>

// Constructor.
Frame::Frame( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
// Destructor.
Frame::~Frame() {}

int Frame::resample(
    CKD *ckd, // Calibration key data.
    int fct // Enahncement factor.
)
{
    size_t dim_spec_enhance = fct*(dim_spec_truth-1)+1;
    vector<double> intens_resample(ckd->dim_fov*dim_spec_enhance);

    // Interpolate onto CKD FOVs. Both the source and destination
    // abscissae are equidistant and start with zero and end with the
    // same number.
    size_t ivp = 0;
    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        // Floating-point index.
        double idx =
          static_cast<double>(ifov) /
          (static_cast<double>(ckd->fov_nfov_vp[ivp]) - 1.0) *
          (static_cast<double>(dim_spat_truth) - 1.0)
          ; // This should move equidistantly from zero to N-1 where N
            // is the input domain size. There are nfov samples.
        size_t intidx = floor(idx);
        if (intidx == dim_spat_truth-1) intidx--; // So that regular interpolation can take place.
        double weightright = idx - intidx;
        double weightleft = 1.0 - weightright;
        size_t &leftidx = intidx;
        size_t rightidx = leftidx+1;

        // If the input is not enhanced, this is not the most efficient
        // implementation in the world, involving an interpolation to
        // exactly the same target. But otherwise, I am juggling with
        // pointers and hopefully, this is still a very small part of
        // the calculation time.
        vector<double> intermediate(dim_spec_truth);
        for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
            intermediate[ispec] =
              weightleft * intens[ivp * dim_spat_truth * dim_spec_truth
                                  + leftidx * dim_spec_truth + ispec]
              + weightright * intens[
                ivp * dim_spat_truth * dim_spec_truth
                + rightidx * dim_spec_truth + ispec];
        }
        linear_interpol(
          dim_spec_truth, // Size of the original array.
          dim_spec_enhance, // Size of the array of interpolated values.
          wavelength.data(), // Independent variable of original array.
          wavelength.data(), // Independent variable of target array.
          intermediate.data(), // Dependent variable of original array.
          &intens_resample[ifov*dim_spec_enhance]); // Dependent variable of target array (output).
    }

    // Overwrite truth with new truth.
    dim_spec_truth = dim_spec_enhance;
    intens = intens_resample;
    return 0;
}

int Frame::convert_units( // {{{
    CKD *ckd
)
{

    // Source unit: photons / (m^2 nm s sr)
    // Destination unit: W / (m^2 um sr)
    // Conversion is 1000 times the energy of a photon in Joules.
    const double c = 2.99792458e17; // Nanometers per second.
    const double h = 6.62607004e-34; // Joule seconds.
    vector<double> conv(dim_spec_truth);
    for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
        conv[ispec] = 1.0e3*h*c/wavelength[ispec];
    }
    double *intens_cur = intens.data();
    size_t ivp = 0;
    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        if (!ckd->vp_mask[ivp]) {
            for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
                intens_cur[ispec] *= conv[ispec];
            }
        }
        intens_cur += dim_spec_truth;
    }

    return 0;

} // }}}

int Frame::read_l1x( // {{{
    size_t iframe,
    CKD *ckd,
    L1X_inputfile *l1x_inputfile
)
{

    l1x_iframe = l1x_inputfile->subset[iframe];

    NcGroup grpnd;
    vector<size_t> strt = {l1x_iframe};
    vector<size_t> cnt = {1};

    l1x_t il1x_start = l1x_inputfile->il1x_start;

    // True scene.
    if (il1x_start == L1X_TRUTH) {
        // Copy fixed cotents: Dimensions and wavelengths.
        dim_spat_truth = l1x_inputfile->dim_spat_truth;
        dim_spec_truth = l1x_inputfile->dim_spec_truth;
        wavelength = l1x_inputfile->wavelength; // Copy the contents.
        // Read raw spectra.
        vector<size_t> strt = { l1x_iframe, 0, 0 };
        vector<size_t> cnt = { 1, l1x_inputfile->dim_spat_truth, l1x_inputfile->dim_spec_truth };
        intens.resize(ckd->dim_vp*l1x_inputfile->dim_spat_truth*l1x_inputfile->dim_spec_truth);
        vector<double *> targets = { intens.data() };
        NcVar var_spectra;
        var_spectra = l1x_inputfile->nc->ncid->getVar("radiance");
        netcdf_check(l1x_inputfile->nc,var_spectra.getVar(strt, cnt, intens.data()));
    }
    return 0;
}

int Frame::write_truth_l1b( // {{{
    CKD *ckd,
    Settings_l1b *set
)
{
    // Write time.
    {
        vector<size_t> strt = {l1x_iframe};
        vector<size_t> cnt = {1};
        netcdf_check(set->nc_l1b,set->var_image_time.putVar(strt,cnt,&image_time));
    }

    // This routine convolves the Gauss with the truth over the input
    // wavelength domain. That is only one wavelength domain for all
    // spectra.

    // Wavelength differences.
    vector<double> dw(dim_spec_truth);
    // In principle, the delta wavelength is defined as half the
    // distance between left and right neighbours. Only at the sides,
    // you miss one of them, and then, we just take whole the distance
    // to the existing neighbour.
    for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
        if (ispec == 0) {
            dw[ispec] = wavelength[ispec+1] - wavelength[ispec];
        } else if (ispec == dim_spec_truth-1) {
            dw[ispec] = wavelength[ispec] - wavelength[ispec-1];
        } else {
            dw[ispec] = 0.5 * (wavelength[ispec+1] - wavelength[ispec-1]);
        }
    }

    // Calculate the Gausses and apply the delta wavelength.
    vector<double> gausses(set->dim_pol_wave*dim_spec_truth);
    double *gauss_cur = gausses.data();
    for (size_t iwave=0 ; iwave<set->dim_pol_wave ; iwave++) {
        double fwhm = set->l1b_wavelength[iwave] / set->resolving_power;
        double sig = fwhm / (2.0*sqrt(2.0*log(2.0)));
        double norm = 0.0;
        for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
            if (abs(wavelength[ispec]-set->l1b_wavelength[iwave])/fwhm > set->gauss_range) gauss_cur[ispec] = 0.0;
            else {
                gauss_cur[ispec] = exp(-pow(wavelength[ispec]-set->l1b_wavelength[iwave],2.0) / (2.0*pow(sig,2.0))) / (sig*sqrt(2.0*PI)) * dw[ispec];
                norm += gauss_cur[ispec];
            }
        }
        writelog(log_debug,"Gauss for wavelength %.7f has normalization %.7f.",set->l1b_wavelength[iwave],norm);
        for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) gauss_cur[ispec] /= norm;
        gauss_cur += dim_spec_truth;
    }

    double *intens_cur = intens.data();
    size_t ivp = 0;
    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        vector<double> i(set->dim_int_wave);
        vector<double> i_polsample(set->dim_pol_wave,0.0);
        vector<double> q(set->dim_pol_wave,0.0);
        vector<double> u(set->dim_pol_wave,0.0);
        vector<double> dolp(set->dim_pol_wave,0.0);
        vector<double> aolp(set->dim_pol_wave,0.0);
        if (ckd->vp_mask[ivp]) {
            for (size_t iwave=0 ; iwave<set->dim_int_wave ; iwave++) i[iwave] = NC_FILL_DOUBLE;
            for (size_t iwave=0 ; iwave<set->dim_pol_wave ; iwave++) {
                i_polsample[iwave] = NC_FILL_DOUBLE;
                q[iwave] = NC_FILL_DOUBLE;
                u[iwave] = NC_FILL_DOUBLE;
                dolp[iwave] = NC_FILL_DOUBLE;
                aolp[iwave] = NC_FILL_DOUBLE;
            }
        } else {
            // True radiance is just an interpolation. That is not even
            // FOV-dependent.
            linear_interpol(
                dim_spec_truth, // Size of the original array.
                set->dim_int_wave, // Size of the array of interpolated values.
                wavelength.data(), // Independent variable of original array.
                set->intensity_wavelength.data(), // Independent variable of target array.
                intens_cur, // Dependent variable of original array.
                i.data() // Dependent variable of target array (output).
            );
            // Now, do all the Gauss convolutions.
            double *gauss_cur = gausses.data();
            for (size_t iwave=0 ; iwave<set->dim_pol_wave ; iwave++) {
                for (size_t ispec=0 ; ispec<dim_spec_truth ; ispec++) {
                    i_polsample[iwave] += intens_cur[ispec] * gauss_cur[ispec];
                }
                gauss_cur += dim_spec_truth;
            }
        }

        // Write them to the NetCDF.
        vector<size_t> strt_int = {l1x_iframe,ifov,0};
        vector<size_t> cnt_int = {1,1,set->dim_int_wave};
        netcdf_check(set->nc_l1b,set->var_radiance_raw.putVar(strt_int,cnt_int,i.data()));
        vector<size_t> strt_pol = {l1x_iframe,ifov,0};
        vector<size_t> cnt_pol = {1,1,set->dim_pol_wave};
        netcdf_check(set->nc_l1b,set->var_intens.putVar(strt_pol,cnt_pol,i_polsample.data()));
        // Progress pointers to next FOV.
        intens_cur += dim_spec_truth;

    }

    return 0;

} // }}}

int Frame::isrf_interpolate(CKD *ckd, Settings_isrf *set)
{
    // Convolve the intensities with an ISRF from the line-by-line
    // grid onto the target grids taken from the CKD.
    vector<double> intens_interp(ckd->dim_fov * ckd->dim_detector_spec, 0.0);
    const double sigma { set->fwhm_gauss / (2.0 * sqrt(2.0 * log(2.0))) };
    const double sigma_inv { 1.0 / (2.0 * sigma * sigma) } ;
    const double norm_inv { 1.0 / (sigma * std::sqrt(2 * std::numbers::pi)) };
    for (int i_fov {}; i_fov < ckd->dim_fov; ++i_fov) {
        // Outer and inner loops are over the target and source wavelength grids
        for (int i_target { i_fov * static_cast<int>(ckd->dim_detector_spec) };
             i_target < (i_fov + 1) * static_cast<int>(ckd->dim_detector_spec);
             ++i_target) {
            for (int i_src {}; i_src < dim_spec_truth; ++i_src) {
                double conv { ckd->wave_target[i_target] - wavelength[i_src] };
                // Reduce the computational cost by considering the
                // limited extent of the Gaussian.
                constexpr double wave_rel_threshold { 3.0 };
                if (std::abs(conv) > wave_rel_threshold * set->fwhm_gauss) {
                    continue;
                }
                conv *= conv * sigma_inv;
                conv = norm_inv * std::exp(-conv);
                double delta_lambda;
                if (i_src == 0) {
                    delta_lambda =
                      wavelength[i_src + 1] - wavelength[i_src];
                } else if (i_src == dim_spec_truth - 1) {
                    delta_lambda = wavelength[dim_spec_truth - 1]
                                   - wavelength[dim_spec_truth - 2];
                } else {
                    delta_lambda =
                      0.5 * (wavelength[i_src + 1] - wavelength[i_src - 1]);
                }
                conv *= delta_lambda;
                intens_interp[i_target] +=
                  conv * intens[i_fov * dim_spec_truth + i_src];
            }
        }
    }
    intens = intens_interp;
    dim_spec_truth = ckd->dim_detector_spec;
    return 0;
}

int Frame::interpolate_truth( // {{{
    CKD *ckd
)
{
    // This routine interpolates the true quantity intens to the CKD
    // wavelengths. This uses the wave_target from the CKD. The member
    // variables will be replaced with the interpolated ones.
    vector<double> intens_interpolated(ckd->dim_fov*ckd->dim_detector_spec);
    // Running pointers.
    vector<double *> sources = {
        intens.data()
    };
    vector<double *> targets = {
        intens_interpolated.data()
    };
    double *wave_target_cur = ckd->wave_target.data();

    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        linear_interpol(
          dim_spec_truth, // Size of the original array.
          ckd->dim_detector_spec, // Size of the array of interpolated values.
          wavelength.data(), // Independent variable of original array.
          wave_target_cur, // Independent variable of target array.
          sources[0], // Dependent variable of original array.
          targets[0]); // Dependent variable of target array (output).
        sources[0] += dim_spec_truth;
        targets[0] += ckd->dim_detector_spec;
        wave_target_cur += ckd->dim_detector_spec;
    }
    intens = intens_interpolated;
    return 0;

} // }}}

int Frame::uncalibrate_spectra( // {{{
    CKD *ckd
)
{
    // Unapply radiometric correction
    double *rad_spectra_cur = ckd->rad_spectra.data();
    double *intens_cur = intens.data();
    size_t ivp = 0;
    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        if (!ckd->vp_mask[ivp]) {
            for (size_t iel=0 ; iel<ckd->dim_detector_spec ; iel++) {
                intens_cur[iel] *= exposure_time / rad_spectra_cur[iel];
            }
        }
        rad_spectra_cur += ckd->dim_detector_spec;
        intens_cur += ckd->dim_detector_spec;
    }
    return 0;

} // }}}

int Frame::draw_on_detector( // {{{
    CKD *ckd
)
{
    image.resize(ckd->npix, 0.0);
    constexpr bool use_splines { true };
    std::vector<double> x_values(ckd->dim_fov);
    std::vector<double> y_values(ckd->dim_fov);
    // Generate an approximate smooth image using splines
    for (int i_spec {}; i_spec < ckd->dim_detector_spec; ++i_spec) {
        for (int i_fov {}; i_fov < ckd->dim_fov; ++i_fov) {
            const int fov_idx { static_cast<int>(ckd->dim_fov - 1
                                                 - i_fov) };
            x_values[i_fov] =
              ckd->fov_ispat[fov_idx * ckd->dim_detector_spec + i_spec];
            y_values[i_fov] =
              intens[fov_idx * ckd->dim_detector_spec + i_spec];

        }
        // LinearSpline spline { x_values, y_values };
        CubicSpline spline { x_values, y_values };
        for (int i_spat {}; i_spat < ckd->dim_detector_spat; ++i_spat) {
            image[i_spat * ckd->dim_detector_spec + i_spec] =
              spline.eval(i_spat);
        }
    }
    // Generate more accurate values for all pixels
    for (int ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
        for (int ifov=0 ; ifov < ckd->dim_fov; ++ifov) {
            size_t *fov_ipix1_cur =
              &ckd->fov_ipix1[ckd->fov_iel_start[ifov]];
            size_t *fov_ipix2_cur =
              &ckd->fov_ipix2[ckd->fov_iel_start[ifov]];
            double *fov_weight1_cur =
              &ckd->fov_weight1[ckd->fov_iel_start[ifov]];
            const int ipix_left = fov_ipix1_cur[ispec];
            const int ipix_right = fov_ipix2_cur[ispec];
            double &weightleft = fov_weight1_cur[ispec];
            const auto s { intens[ifov*ckd->dim_detector_spec + ispec] };
            if (std::abs(image[ipix_left]) < 1e-100) {
                image[ipix_left] = s;
            }
            image[ipix_right] = (s - weightleft * image[ipix_left]) / (1.0 - weightleft);
        }
    }
    return 0;
} // }}}

int Frame::apply_straylight(const Settings_main& set, CKD *ckd)
{
    for (int i {}; i < static_cast<int>(image.size()); ++i) {
        if (std::isnan(image[i])) {
            image[i] = image[i - 1];
        }
    }
    std::vector<double> conv_result(ckd->npix, 0.0);
    for (int i_kernel {}; i_kernel < ckd->stray.n_kernels; ++i_kernel) {
        std::vector<double> image_weighted { image };
        for (int i {}; i < ckd->npix; ++i) {
            image_weighted[i] *= ckd->stray.weights[i_kernel][i];
        }
        const int image_n_rows {
            ckd->stray.edges[i_kernel * box::n + box::t]
            - ckd->stray.edges[i_kernel * box::n + box::b]
        };
        const int image_n_cols {
            ckd->stray.edges[i_kernel * box::n + box::r]
            - ckd->stray.edges[i_kernel * box::n + box::l]
        };
        std::vector<double> sub_image(image_n_rows * image_n_cols);
        for (int i {}; i < image_n_rows; ++i) {
            for (int j {}; j < image_n_cols; ++j) {
                sub_image[i * image_n_cols + j] =
                  image_weighted
                  [(i + ckd->stray.edges[i_kernel * box::n + box::b])
                   * ckd->dim_detector_spec
                   + j
                   + ckd->stray.edges[i_kernel * box::n + box::l]];
            }
        }
        std::vector<double> conv_result_sub {};
        convolve_fft(image_n_rows,
                     image_n_cols,
                     sub_image,
                     ckd->stray.kernel_rows[i_kernel],
                     ckd->stray.kernel_cols[i_kernel],
                     ckd->stray.kernel_fft_sizes[i_kernel],
                     ckd->stray.kernels_fft[i_kernel],
                     conv_result_sub);
        for (int i {}; i < image_n_rows; ++i) {
            for (int j {}; j < image_n_cols; ++j) {
                conv_result
                  [(i + ckd->stray.edges[i_kernel * box::n + box::b])
                   * ckd->dim_detector_spec + j
                   + ckd->stray.edges[i_kernel * box::n + box::l]] +=
                  conv_result_sub[i * image_n_cols + j];
            }
        }
    }
    std::vector<double> image_conv(ckd->npix);
    for (int i {}; i < ckd->npix; ++i) {
        image_conv[i] =
          (1.0 - ckd->stray.eta[i]) * image[i] + conv_result[i];
    }
    std::swap(image_conv, image);
    return 0;
}

int Frame::apply_prnu(CKD *ckd)
{
    // Pixel-response non-uniformity
    for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
        if (!ckd->mask[ipix]) {
            image[ipix] *= ckd->prnu_prnu[ipix];
        }
    }
    return 0;
}

int Frame::apply_nonlinearity(CKD *ckd)
{
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
    double *nonlin_fit_cur = &ckd->nonlin_fit[iexptime_best*ckd->dim_nonlin_spline];
    double *nonlin_signal_scale_cur = &ckd->nonlin_signal_scale[iexptime_best];
    for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
        if (!ckd->mask[ipix]) {
            double abscissa = image[ipix] / (*nonlin_signal_scale_cur);
            handle(b.jacapply(1,&abscissa,nonlin_fit_cur,&image[ipix]));
        }
        nonlin_fit_cur += ckd->dim_nonlin_spline*ckd->dim_nonlin_exptime;
        nonlin_signal_scale_cur += ckd->dim_nonlin_exptime;
    }
    return 0;
}

int Frame::apply_dark_current(CKD *ckd)
{

    // This one only writes the image with dark current.
    image_with_current = image; // Copy contents.

    // Evaluate the fit. The B-splines are fixed to a polynomial where the first term
    // is the value at nominal temperature.
    vector<double> knots = {ckd->dark_nominal_temperature,ckd->dark_nominal_temperature+1.0};
    Bspline b(this,ckd->dim_dark_order,2,knots.data());
    vector<double> terms(ckd->dim_dark_order);

    handle(b.jaccalc(1,&detector_temperature,terms.data()));
    for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
        if (!ckd->mask[ipix]) {
            for (size_t iorder=0 ; iorder<ckd->dim_dark_order ; iorder++) {
                image_with_current[ipix] += ckd->dark_current[ipix*ckd->dim_dark_order+iorder]*terms[iorder] * exposure_time; // Add dark current.
            }
        }
    }
    return 0;
}

int Frame::apply_dark_offset(CKD *ckd)
{
    // Evaluate the fit. The B-splines are fixed to a polynomial where the first term
    // is the value at nominal temperature.
    vector<double> knots = {ckd->dark_nominal_temperature,ckd->dark_nominal_temperature+1.0};
    Bspline b(this,ckd->dim_dark_order,2,knots.data());
    vector<double> terms(ckd->dim_dark_order);

    handle(b.jaccalc(1,&detector_temperature,terms.data()));
    image = image_with_current;
    for (size_t ipix=0 ; ipix<ckd->npix ; ipix++) {
        if (!ckd->mask[ipix]) {
            for (size_t iorder=0 ; iorder<ckd->dim_dark_order ; iorder++) {
                image[ipix] += ckd->dark_offset[ipix*ckd->dim_dark_order+iorder]*terms[iorder]; // Add dark offset.
            }
        }
    }
    return 0;
}
