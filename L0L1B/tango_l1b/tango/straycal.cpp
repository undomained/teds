// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "ckd.h"
#include "l1a.h"
#include "straycal.h"

#include <algorithm>
#include <array>
#include <complex>
#include <fftw3.h>
#include <cassert>

// Settings functions.
Settings_straycal::Settings_straycal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "stray";
} // }}}
Settings_straycal::~Settings_straycal() {}
int Settings_straycal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(padding);
    recognize_setting(moving_kernel_intensity_threshold);
    recognize_setting(dry_run);

    return 0;

} // }}}

static void rotateKernel( // {{{
    const int kernel_n_rows,
    const int kernel_n_cols,
    vector<double> &kernel
)
{
    const int kernel_n_rows_half { kernel_n_rows/2 };
    const int kernel_n_cols_half { kernel_n_cols/2 };
    const int kernel_n_rows_half_one { kernel_n_rows_half + 1 };
    const int kernel_n_cols_half_one { kernel_n_cols_half + 1 };
    vector<double> kernel_initial(kernel);
    for (int i {}; i < kernel_n_rows_half_one; ++i) {
        const int idx_i { (i+kernel_n_rows_half)*kernel_n_cols };
        const int idx_f { i*kernel_n_cols };
        for (int j {}; j < kernel_n_cols_half_one; ++j) {
            kernel[idx_f+j] = kernel_initial[idx_i+j+kernel_n_cols_half];
        }
    }
    for (int i { kernel_n_rows_half_one }; i < kernel_n_rows; ++i) {
        const int idx_i { (i-kernel_n_rows_half_one)*kernel_n_cols };
        const int idx_f { i*kernel_n_cols };
        for (int j { kernel_n_cols_half_one }; j < kernel_n_cols; ++j) {
            kernel[idx_f+j] = kernel_initial[idx_i+j-kernel_n_cols_half_one];
        }
    }
    for (int i { kernel_n_rows_half_one }; i < kernel_n_rows; ++i) {
        const int idx_i { (i-kernel_n_rows_half_one)*kernel_n_cols };
        const int idx_f { i*kernel_n_cols };
        for (int j {}; j < kernel_n_cols_half_one; ++j) {
            kernel[idx_f+j] = kernel_initial[idx_i + j+kernel_n_cols_half];
        }
    }
    for (int i {}; i < kernel_n_rows_half_one; ++i) {
        const int idx_i { (i+kernel_n_rows_half)*kernel_n_cols };
        const int idx_f { i*kernel_n_cols };
        for (int j { kernel_n_cols_half_one }; j < kernel_n_cols; ++j) {
            kernel[idx_f+j] = kernel_initial[idx_i + j-kernel_n_cols_half_one];
        }
    }
} // }}}

// Constructor for the Straycal CKD structure.
Straycal::Straycal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("stray");
    set = make_unique<Settings_straycal>(this);
    Processor::set = set.get();
    own_skip = &ckd->stray_skip;
} // }}}
Straycal::~Straycal() {}

// Straylight correction calibration software.
int Straycal::process_batch( // {{{
    size_t ibatch // Meaningless zero here.
)
{
    // PART 0 - Dimensions.
    static constexpr double default_relative_padding { 0.1 };
    ckd->stray_kernel_n_rows = 2*ckd->dim_detector_spat + 1;
    ckd->stray_kernel_n_cols = 2*ckd->dim_detector_spec + 1;
    const int kernel_n_rows { static_cast<int>(ckd->stray_kernel_n_rows) };
    const int kernel_n_cols { static_cast<int>(ckd->stray_kernel_n_cols) };
    const int stray_kernel_size { kernel_n_rows*kernel_n_cols };
    // Since we are dealing with a two-dimensional detector, several
    // quantities here are computed over two directions.
    static constexpr int n_dirs { 2 };
    const array<int, n_dirs> kernel_mid { kernel_n_rows/2, kernel_n_cols/2 };

    // PART 1 - Shift all calibration images to their centers
    vector<vector<double>> I_centered(nl1a);
    vector<array<int, n_dirs>> centers(nl1a);
    vector<double> I_cal(ckd->dim_detector_spat*ckd->dim_detector_spec);
    for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
        L1A *l1a = l1a_instances[i_cal];
        copy(l1a->image, l1a->image + I_cal.size(), I_cal.begin());
        for (int i {}; i < static_cast<int>(I_cal.size()); ++i) {
            if (l1a->pixelmask[i]) I_cal[i] = 0.0;
        }
        I_centered[i_cal] = vector<double>(stray_kernel_size, 0.0);
        array<int, n_dirs> max_coords {};
        double max_value {};
        for (int m {}; m < static_cast<int>(ckd->dim_detector_spat); ++m) {
            for (int n {}; n < static_cast<int>(ckd->dim_detector_spec); ++n) {
                if (I_cal[m*ckd->dim_detector_spec+n] > max_value) {
                    max_value = I_cal[m*ckd->dim_detector_spec+n];
                    max_coords = { m, n };
                }
            }
        }
        copy(max_coords.begin(), max_coords.end(), centers[i_cal].begin());
        for (int m {}; m < static_cast<int>(ckd->dim_detector_spat); ++m) {
            const int m_center { m+kernel_mid[0]-max_coords[0] };
            const int mc_idx { m_center*kernel_n_cols };
            for (int n {}; n < static_cast<int>(ckd->dim_detector_spec); ++n) {
                const int n_center { n+kernel_mid[1]-max_coords[1] };
                const int n_cols { static_cast<int>(ckd->dim_detector_spec) };
                const int I_cal_idx { m*n_cols+n };
                I_centered[i_cal][mc_idx+n_center] = I_cal[I_cal_idx];
            }
        }
    }

    // PART 2 - Create a constant kernel using the lowest values. Here
    //          ckd->stray_kernel is actually the full kernel. It
    //          becomes the constant part of the stray light kernel in
    //          PART 6.
    vector<double> constant_kernel(stray_kernel_size, 0.0);
    vector<double> values_at_pixel(nl1a);
    for (int m {}; m < kernel_n_rows; ++m) {
        for (int n {}; n < kernel_n_cols; ++n) {
            // Use the lowest value at point that is greater than zero
            for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
                const int c_idx { m*kernel_n_cols+n };
                values_at_pixel[i_cal] = I_centered[i_cal][c_idx];
            }
            sort(values_at_pixel.begin(), values_at_pixel.end());
            // Find the lowest positive value
            constant_kernel[m*ckd->stray_kernel_n_cols+n] = values_at_pixel[0];
        }
    }

    // PART 3 - Capture all moving kernels by subtracting the full
    //          kernel from each calibration measurement. In
    //          principle, what is left is equal to the moving
    //          kernel. Also capture the positions of the moving
    //          kernel centers.
    vector<vector<double>> K_moves(nl1a);
    vector<array<double, n_dirs>> K_centers(nl1a);
    for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
        K_moves[i_cal] = vector<double>(stray_kernel_size, 0.0);
        for (int i {}; i < static_cast<int>(constant_kernel.size()); ++i) {
            K_moves[i_cal][i] = I_centered[i_cal][i] - constant_kernel[i];
        }
        // After subtracting the full kernel, everything that is
        // within moving_kernel_intensity_threshold of the maximum
        // value defines a moving kernel.
        double max_value {};
        for (const auto &el : K_moves[i_cal]) max_value = max(max_value, el);
        double intensity_sum {};
        K_centers[i_cal] = {};
        for (int m {}; m < kernel_n_rows; ++m) {
            for (int n {}; n < kernel_n_cols; ++n) {
                if (K_moves[i_cal][m*kernel_n_cols+n]
                    >= set->moving_kernel_intensity_threshold*max_value) {
                    intensity_sum += K_moves[i_cal][m*kernel_n_cols+n];
                    // Moving kernel centers are determined by the
                    // center of mass formula. If the intensity is
                    // uniform, it is the same as averaging over all
                    // pixel positions that make up the moving kernel.
                    K_centers[i_cal][0] += m*K_moves[i_cal][m*kernel_n_cols+n];
                    K_centers[i_cal][1] += n*K_moves[i_cal][m*kernel_n_cols+n];
                }
            }
        }
        assert(intensity_sum > 0.0 && "Moving kernel has zero intensity");
        for (int i_dir {}; i_dir < n_dirs; ++i_dir) {
            K_centers[i_cal][i_dir] = K_centers[i_cal][i_dir]/intensity_sum;
            K_centers[i_cal][i_dir] -= kernel_mid[i_dir];
        }
    }

    // PART 4 - Find the alpha and beta parameters using linear
    //          regression. Moving kernel movement is described by
    //          y_i(x) = beta_i - alpha_i*x_i,
    //          where i is the direction (across rows or columns).
    array<double, n_dirs> alpha_means {};
    array<double, n_dirs> beta_means {};
    for (int i_dir {}; i_dir < n_dirs; ++i_dir) {
        for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
            alpha_means[i_dir] += centers[i_cal][i_dir];
            beta_means[i_dir] += K_centers[i_cal][i_dir];
        }
    }
    for (auto &el : alpha_means) el /= nl1a;
    for (auto &el : beta_means) el /= nl1a;
    // Find alphas
    //
    //         Sum_i (x_i - x_mean)(y_i - y_mean)
    // alpha = ----------------------------------,
    //               Sum_i (x_i - x_mean)^2
    //
    // where x_i are centers of the calibration signal (position of
    // point light source on the detector) and y_i are centers of the
    // moving kernels.
    array<double, n_dirs> alphas {};
    array<double, n_dirs> squares {};  // (x_i - x_mean)^2
    for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
        for (int i_dir {}; i_dir < n_dirs; ++i_dir) {
            const double diff { centers[i_cal][i_dir] - alpha_means[i_dir] };
            alphas[i_dir] += diff*(K_centers[i_cal][i_dir] - beta_means[i_dir]);
            squares[i_dir] += diff*diff;
        }
    }
    for (int i_dir {}; i_dir < n_dirs; ++i_dir) alphas[i_dir] /= squares[i_dir];
    // Find betas
    //
    // beta = y_mean - alpha*x_mean
    array<double, n_dirs> betas {};
    for (int i_dir {}; i_dir < n_dirs; ++i_dir) {
        betas[i_dir] = -(beta_means[i_dir] - alphas[i_dir]*alpha_means[i_dir]);
    }

    // PART 5 - Construct the final moving kernel by moving all moving
    //          kernels to the centers and taking the average.
    const auto roundToInt { [](double x) { return static_cast<int>(round(x)); } };
    vector<double> moving_kernel(stray_kernel_size, 0.0);
    for (int i_cal {}; i_cal < static_cast<int>(nl1a); ++i_cal) {
        for (int m {}; m < kernel_n_rows; ++m) {
            const int center_m { m + roundToInt(K_centers[i_cal][0]) };
            if (center_m >= 0 && center_m < kernel_n_rows) {
                for (int n {}; n < kernel_n_cols; ++n) {
                    const int center_n { n + roundToInt(K_centers[i_cal][1]) };
                    if (center_n >= 0 && center_n < kernel_n_cols) {
                        const int idx { center_m*kernel_n_cols + center_n };
                        moving_kernel[m*kernel_n_cols+n] += K_moves[i_cal][idx];
                    }
                }
            }
        }
    }
    for (auto &el : moving_kernel) el /= nl1a;

    // PART 6 - Normalize kernels and create the constant kernel by
    //          removing the central element of the total kernel.
    double K_total_sum {};
    for (const auto &el : constant_kernel) K_total_sum += el;
    double K_move_centered_sum {};
    for (const auto &el : moving_kernel) K_move_centered_sum += el;
    const double total_intensity { K_total_sum + K_move_centered_sum };
    for (auto &el : constant_kernel) el /= total_intensity;
    for (auto &el : moving_kernel) el /= total_intensity;
    K_move_centered_sum = 0.0;
    for (const auto &el : moving_kernel) K_move_centered_sum += el;
    // Now constant_kernel becomes the constant kernel
    constant_kernel[kernel_mid[0]*kernel_n_cols + kernel_mid[1]] = 0.0;
    // Internal scattering factor
    ckd->stray_eta = 0.0;
    for (const auto &el : constant_kernel) ckd->stray_eta += el;
    ckd->stray_eta += K_move_centered_sum;

    // PART 7 - Rotate the kernels and store their Fourier transforms.
    rotateKernel (
        ckd->stray_kernel_n_rows,
        ckd->stray_kernel_n_cols,
        constant_kernel
    );
    rotateKernel (
        ckd->stray_kernel_n_rows,
        ckd->stray_kernel_n_cols,
        moving_kernel
    );
    // Note that because the matrices are real, only half the elements
    // need to be stored.
    ckd->stray_kernel_fft_size = kernel_n_rows*(kernel_n_cols/2+1);
    // First do the constant part of the kernel
    ckd->stray_kernel_fft.resize(ckd->stray_kernel_fft_size);
    fftw_plan fft_plan = fftw_plan_dft_r2c_2d (
        ckd->stray_kernel_n_rows,
        ckd->stray_kernel_n_cols,
        constant_kernel.data(),
        reinterpret_cast<fftw_complex*>(ckd->stray_kernel_fft.data()),
        FFTW_ESTIMATE
    );
    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);
    // Now the moving kernel
    ckd->stray_moving_kernel_fft.resize(ckd->stray_kernel_fft_size);
    fft_plan = fftw_plan_dft_r2c_2d (
        ckd->stray_kernel_n_rows,
        ckd->stray_kernel_n_cols,
        moving_kernel.data(),
        reinterpret_cast<fftw_complex*>(ckd->stray_moving_kernel_fft.data()),
        FFTW_ESTIMATE
    );
    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);

    // PART 8 - Create the auxiliary array for transforming
    //          matrices. These are used in L1A::calibrate_detector
    //          when computing the convolution with a moving kernel.
    if (set->padding == -1) {
        set->padding = default_relative_padding*ckd->dim_detector_spat;
    }
    ckd->stray_transformed_n_rows = ckd->dim_detector_spat + 2*set->padding;
    ckd->stray_transformed_n_cols = ckd->dim_detector_spec + 2*set->padding;
    const int tr_n_rows { static_cast<int>(ckd->stray_transformed_n_rows) };
    const int tr_n_cols { static_cast<int>(ckd->stray_transformed_n_cols) };
    const int stray_transformed_size { tr_n_rows*tr_n_cols };
    ckd->stray_transform_indices.resize(stray_transformed_size);
    ckd->stray_transform_deltas.resize(4*stray_transformed_size);
    auto &transform_indices { ckd->stray_transform_indices };
    auto &transform_deltas { ckd->stray_transform_deltas };
    const double norm { (1+alphas[0])*(1+alphas[1]) };
    int delta_idx { -1 };
    for (int i { -set->padding }; i < tr_n_rows-set->padding; ++i) {
        const double i_new_d { (i+betas[0])/(1+alphas[0]) };
        const int i_new { static_cast<int>(i_new_d) };
        const double i_delta { i_new_d - i_new };
        for (int j { -set->padding }; j < tr_n_cols-set->padding; ++j) {
            const double j_new_d { (j+betas[1])/(1+alphas[1]) };
            const int j_new { static_cast<int>(j_new_d) };
            const double j_delta { j_new_d - j_new };
            const int idx { (i+set->padding)*tr_n_cols+j+set->padding };
            transform_indices[idx] = i_new*ckd->dim_detector_spec + j_new;
            transform_deltas[++delta_idx] = (1-i_delta)*(1-j_delta)/norm;
            transform_deltas[++delta_idx] = i_delta*(1-j_delta)/norm;
            transform_deltas[++delta_idx] = (1-i_delta)*j_delta/norm;
            transform_deltas[++delta_idx] = i_delta*j_delta/norm;
        }
    }

    // PART 9 - If dry_run is true, still do all the computations but
    //          don't apply the stray light correction in l1a.
    ckd->stray_dry_run = set->dry_run;

    return 0;
} // }}}
