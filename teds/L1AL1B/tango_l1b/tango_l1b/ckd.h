// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "logger.h"
#include "settings_proc.h" // For members with class Calibration_options.

#include <complex>

namespace tango {

// Forward declaration.
class NetCDF_object;
class Settings_main;

class CKD : public Logger { // {{{

    // All data is public. Think of it like a struct.
    public:
    // Constructor.
    CKD(
        Logger *creator
    );
    ~CKD();

    // Level of the process for which this CKD is mature enough for.
    level_t lev; // To be set during read routine.

    // Output NetCDF object.
    bool write_ckd; // Flag for writing an output CKD file.
    unique_ptr<NetCDF_object> nc_ckd; // To prevent the .get() expressions all the time.

    // Dimennsion sizes.
    // Level DIM. These dimensions should always be known. The first step
    // is to acquire these dimensions. Then they are known at any moment
    // something serious happens.
    size_t dim_detector_spec; // Number of detector pixels in the spectral direction.
    size_t dim_detector_spat; // Number of detector pixels in the spatial direction.
    size_t dim_vp; // Number of view ports.
    // Level DARK. Dark calibration temperature dependence fitting order.
    size_t dim_dark_order; // This is the dimension, the same as the B-spline order. A second order polynomial is B-spline order 3, and also dim_dark_order is 3.
    // Level NONLIN. Non-linearity correction.
    size_t dim_nonlin_exptime; // Number of exposure times for which non-linearity is calculated.
    size_t dim_nonlin_spline; // Number of B-splines in non-lineary correction fit.
    size_t dim_nonlin_knot; // Number of knots in the B-spline of the non-linearity fit.
    // Fields of view, produced by fovcal, but relevant for later steps as well.
    size_t dim_fov = 0; // Number of spectra in one swath. Default needed because of partial L1B processing.

    // Derived dimension size.
    size_t npix; // Number of pixels in one image, the product of dim_detector_spat and dim_detector_spec.

    // NetCDF dimension identifiers that need to be saved because they are used by
    // different steps.
    NcDim dimid_pol; // This is a fixed-sized one, always two.
    NcDim dimid_vec; // This is a fixed-sized one, always three.
    NcDim dimid_detector_spat;
    NcDim dimid_detector_spec;
    NcDim dimid_vp;
    NcDim dimid_fov;

    // Mask, used by more processors, so saved on structure.
    NcVar var_mask;
    NcVar var_vp_mask;

    // Pixel mask.
    // This pixel mask is updated during all the detector calibration steps.
    // As it is not owned by one particular step, there is no prefix.
    vector<bool> mask; // Dead pixel mask.
    size_t nliving; // Number of living pixels.
    // Viewport mask.
    // This mask can be updated at any step by skipping a viewport. Of
    // course, this makes most sense at a calibration step from FOVcal
    // onwards, because only then, viewports get a meaning.
    // This mask is not used operationally. This is only there to save
    // time in test calculations or test L1B calculations, that you can
    // process or test one viewport without bothering about the others.
    vector<bool> vp_mask; // Skipped viewport mask.

    // Dark CKD.
    bool dark_skip = false; // Flag for skipping dark correction step.
    vector<double> dark_offset; // Fixed dark signal (independent of integration time).
    vector<double> dark_current; // Dark signal added per second of integration time.
    double dark_nominal_temperature; // Nominal temperature temperature-difference fit.
    // Diagnostic datk CKD.
    double *diag_dark_chi2;

    // Noise CKD.
    bool noise_skip = false; // Flag for skipping noise calibration step.
    vector<double> noise_g; // Signal-dependent noise term.
    vector<double> noise_n; // Signal-independent noise term.

    // Nonlin CKD.
    bool nonlin_skip = false; // Flag for skipping non-linearity step.
    size_t nonlin_order; // B-spline order of the non-linearity correction.
    vector<double> nonlin_knots; // B-spline knots for non-linearity correction.
    vector<double> nonlin_fit; // B-spline fit for non-linearity correction to be used.
    vector<double> nonlin_exptimes; // Exposure times for which non-linearity is calculated.
    vector<double> nonlin_signal_scale; // Scaling signal for the abscissa of the non-linearity fit.
    // Diagnostic (useless) nonlin CKD.
    double *diag_nonlin_lin_slope;
    double *diag_nonlin_lin_chi2;
    double *diag_nonlin_chi2;

    // PRNU CKD.
    bool prnu_skip = false; // Flag for skipping PRNU step.
    vector<double> prnu_prnu; // Pixel response non-uniformity.

    // Straylight CKD.
    bool stray_skip = false; // Flag for skipping straylight step.
    bool stray_dry_run;
    size_t stray_kernel_n_rows;
    size_t stray_kernel_n_cols;
    size_t stray_kernel_fft_size;
    vector<complex<double> > stray_kernel_fft;
    vector<complex<double> > stray_moving_kernel_fft;
    size_t stray_transformed_n_rows;
    size_t stray_transformed_n_cols;
    vector<int> stray_transform_indices;
    vector<double> stray_transform_deltas;
    double stray_eta;
    struct
    {
        int n_kernels {};
        // Number of spatial and spectral elements of the unbinned
        // image. This is determined by the detector dimensions and
        // reduced_kernels.
        int n_spatial {};
        int n_spectral {};
        // Number of rows of each kernel
        std::vector<int> kernel_rows {};
        // Number of columns of each kernel
        std::vector<int> kernel_cols {};
        // The FFT array size of each kernel
        std::vector<int> kernel_fft_sizes {};
        // Fourier transforms of the kernels
        std::vector<std::vector<std::complex<double>>> kernels_fft {};
        // Total internal scattering factor
        std::vector<double> eta {};
        // Kernel weights
        std::vector<std::vector<double>> weights {};
        // Boundaries of subimages that must be extracted for the
        // convolutions. The order of coefficients is 'bottom', 'top',
        // 'left', 'right'.
        std::vector<int> edges {};
    } stray;

    // FOV CKD.
    vector<uint32_t> fov_nfov_vp; // Number of fields of view per viewport.
    vector<double> fov_ispat; // Spatial detector pixel coordinate of spectra (only for unbinned CKD).
    vector<size_t> fov_ipix1; // First pixel involved in interpolation for all L1B spectra.
    vector<size_t> fov_ipix2; // Second pixel involved in interpolation for all L1B spectra.
    vector<double> fov_weight1; // Weight factor for the first pixel in this interpolation.
    vector<double> fov_act_angles; // Swath positions for each FOV calibration measurement (In normalized, 0=left, 1=right perspective).
    vector<size_t> fov_iel_start; // Spectral element where a new FOV starts.
    vector<size_t> fov_dims_spec; // Size of the spectra of each FOV.
    size_t dim_fov_spec_total; // Total of the sizes of each FOV.

    // Swath CKD.
    bool swath_skip = false;
    vector<double> swath_swathvectors; // Pointing vector per FOV in satellite coordinates.
    vector<double> swath_vectorplane_normals; // Normals of planes fitted through swath vectors of one viewport.

    // Wave CKD.
    vector<double> wave_spectra; // Spectra of wavelengths calibrated by extracted spectra.
    vector<double> wave_target; // Common wavelength grid for each pair of S+ and S- spectra.

    // Radiometric CKD.
    bool rad_skip = false; // Flag for skipping radiometric calibration.
    vector<double> rad_spectra; // Radiometric calibration factor for extracted spectra.

    // Polarimetric CKD.
    vector<double> pol_m_q; // Mueller matrix elements for Q.
    vector<double> pol_m_u; // Mueller matrix elements for U.
    vector<double> pol_m_t; // Mueller matrix elements for telescope polarization
    double *diag_pol_delta; // Retardance.
    double *diag_pol_eff_a; // Efficiency along one axis (q with zero tilt).
    double *diag_pol_eff_b; // Efficiency along the other axis (u with zero tilt).
    double *diag_pol_tilt; // Tilt angle of the ellipse.

    // Detector options.
    // These are the detector option structures used to generate each of
    // the CKD steps. These are only used to warn for inconsistency and
    // to have them into the CKD file so that the user can do it manually
    // as well.
    // Level dim is totally uncoupled to any interaction with L1A, because
    // without it, there is no detector size or shape, so reading a L1A
    // cannot be done safely (like walking over a plank blindfolded).
    // All other steps have L1A measurements involved, so speaking about
    // detector calibration options makes sense. Possibly, all optional
    // stuff is more advanced than where you are, but then the options
    // are there, but no option is relevant.
    // We will save the forced options as well as the optional options.
    // Of course, writing forced options into the CKD is CKD-file pollution,
    // but it informs the user of what is done and it saves the program
    // the effort for remembering what is forced in what step (and that
    // is error prone when chaning these forced options).
    Calibration_options opt_dark;
    Calibration_options opt_noise;
    Calibration_options opt_nonlin;
    Calibration_options opt_prnu;
    Calibration_options opt_stray;
    Calibration_options opt_fov;
    Calibration_options opt_swath;
    Calibration_options opt_wave;
    Calibration_options opt_rad;
    Calibration_options opt_pol;

    // Member functions.
    public:
    // Reads CKD from open file up to a certain level.
    // For the very first step, nothing is read.
    int read(
        Settings_main &set,
        level_t lev_target,
        bool write
    );

    // Writes the current step into the CKD making the total CKD one step more
    // mature (increasing its level by one).
    // It is assumed that all diagnostic CKD of the current step exists.
    // For L1B and later, this routine does not do anything, it does not even
    // look at the output NetCDF structure.
    int writestep();

    // Checks calibration options.
    int check_opts(
        Calibration_options opt
    );

    private:
    // Reads used calibration option from group attributes. Level is
    // known, but group and calibration option structures are not
    // directly accessible by level number.
    int read_opt(
        NetCDF_object *nc,
        NcGroup grp,
        Calibration_options &opt
    );

    // Writes used calibration optinos to group attributes.
    int write_opt(
        NcGroup grp,
        Calibration_options opt
    );

    // Checks calibration options and saves the one for the actual step.
    int check_opt(
        Calibration_options opt_user,
        level_t lev_user,
        Calibration_options &opt_ref,
        level_t lev_ref,
        string stepname
    );

}; // }}}

} // namespace tango
