#pragma once

#include <tango_l1b/ckd.h>
#include <tango_l1b/fourier.h>
#include <tango_l1b/l1.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <numeric>

// Read an array from the common fixtures folder and return the array
// dimensions (required for constructing the CKD).
template <typename T>
auto readArrayFromTxt(const std::string& filename,
                      std::vector<T>& array) -> std::array<int, 2>
{
    array.clear();
    std::ifstream in { filename };
    std::string line {};
    int n_rows {};
    while (std::getline(in, line)) {
        std::stringstream ss { line };
        T val {};
        while (ss >> val) {
            array.push_back(val);
        }
        ++n_rows;
    }
    return { n_rows, static_cast<int>(array.size()) / n_rows };
}

auto readNonlinCKD(const std::string& fixture_dir,
                   std::vector<double>& x_values,
                   std::vector<double>& y_values) -> void
{
    std::vector<double> buf {};
    readArrayFromTxt(fixture_dir + "/ckd_nonlin.txt", buf);
    x_values.resize(buf.size() / 2);
    y_values.resize(buf.size() / 2);
    for (int i {}; i < static_cast<int>(x_values.size()); ++i) {
        x_values[i] = buf[2 * i];
        y_values[i] = buf[2 * i + 1];
    }
}

auto genStrayCKD(const std::string& fixture_dir,
                 const int n_rows,
                 const int n_cols,
                 tango::CKD& ckd) -> void
{
    ckd.stray.n_kernels = 2;
    std::vector<std::vector<double>> kernels(ckd.stray.n_kernels);
    const auto [n_kernel_rows, n_kernel_cols] { readArrayFromTxt(
      fixture_dir + "/ckd_stray_kernel.txt", kernels.front()) };
    ckd.stray.kernel_rows.resize(ckd.stray.n_kernels);
    ckd.stray.kernel_cols.resize(ckd.stray.n_kernels);
    ckd.stray.kernels_fft.resize(ckd.stray.n_kernels);
    ckd.stray.kernel_fft_sizes.resize(ckd.stray.n_kernels);
    ckd.stray.eta.resize(ckd.stray.n_kernels);
    ckd.stray.weights.resize(ckd.stray.n_kernels);
    ckd.stray.edges.resize(ckd.stray.n_kernels * tango::box::n);
    for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
        kernels[i_kernel] = kernels.front();
        ckd.stray.kernel_rows[i_kernel] = n_kernel_rows;
        ckd.stray.kernel_cols[i_kernel] = n_kernel_cols;
        ckd.stray.kernel_fft_sizes[i_kernel] =
          tango::getFFTSize(n_kernel_rows, n_kernel_cols);
        ckd.stray.kernels_fft[i_kernel].resize(
          ckd.stray.kernel_fft_sizes[i_kernel]);
        tango::fft_r2c(n_kernel_rows,
                       n_kernel_cols,
                       kernels[i_kernel],
                       ckd.stray.kernels_fft[i_kernel]);
        ckd.stray.weights[i_kernel].assign(n_rows * n_cols, 1.0);
        ckd.stray.edges[i_kernel * tango::box::n + tango::box::b] = 0;
        ckd.stray.edges[i_kernel * tango::box::n + tango::box::t] = n_rows;
        ckd.stray.edges[i_kernel * tango::box::n + tango::box::l] = 0;
        ckd.stray.edges[i_kernel * tango::box::n + tango::box::r] = n_cols;
    }
    constexpr double eta { 0.10 };
    ckd.stray.eta.assign(n_rows * n_cols, eta);
}

// Read CKD for tests. Different CKD sections are read one by one from
// text files.
auto readCKD(const std::string& fixture_dir, tango::CKD& ckd) -> void
{
    // Basic dimensions and the pixel mask
    const auto [n_rows, n_cols] { readArrayFromTxt(
      fixture_dir + "/ckd_pixel_mask.txt", ckd.pixel_mask) };
    ckd.n_detector_rows = n_rows;
    ckd.n_detector_cols = n_cols;
    ckd.npix = ckd.n_detector_rows * ckd.n_detector_cols;
    ckd.n_detector_rows_binned = ckd.n_detector_rows;
    ckd.n_detector_cols_binned = ckd.n_detector_cols;
    ckd.npix_binned = ckd.npix;
    // Dark CKD
    readArrayFromTxt(fixture_dir + "/ckd_dark_offset.txt", ckd.dark.offset);
    readArrayFromTxt(fixture_dir + "/ckd_dark_current.txt", ckd.dark.current);
    // Noise CKD
    readArrayFromTxt(fixture_dir + "/ckd_conversion_gain.txt", ckd.noise.g);
    readArrayFromTxt(fixture_dir + "/ckd_read_noise.txt", ckd.noise.n2);
    for (double& val : ckd.noise.n2) {
        val *= val;
    }
    // Nonlinearity CKD
    std::vector<double> x_values {};
    std::vector<double> y_values {};
    readNonlinCKD(fixture_dir, x_values, y_values);
    ckd.nonlin.spline = { x_values, y_values };
    // PRNU CKD
    readArrayFromTxt(fixture_dir + "/ckd_prnu.txt", ckd.prnu.prnu);
    // Stray light CKD
    genStrayCKD(fixture_dir, ckd.n_detector_rows, ckd.n_detector_cols, ckd);
}

auto readL1(const std::string& fixture_dir, tango::L1& l1) -> void
{
    l1.n_alt = 1;
    l1.exposure_time = 0.0460833;
    readArrayFromTxt(fixture_dir + "/detector_image.txt", l1.signal);
    l1.noise.assign(l1.signal.size(), 1.0);
}

auto absSum(std::vector<double> array) -> double
{
    return std::accumulate(
      array.begin(), array.end(), 0.0, [](const double s, const double a) {
          return s + std::abs(a);
      });
}
