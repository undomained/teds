#pragma once

#include <tango_l1b/ckd.h>
#include <tango_l1b/l1.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <fstream>
#include <numeric>

// Read an array from the common fixtures folder and return the array
// dimensions (required for constructing the CKD).
auto readArrayFromTxt(const std::string& filename,
                      std::vector<double>& array) -> std::array<int, 2>
{
    array.clear();
    std::ifstream in { filename };
    std::string line {};
    int n_rows {};
    while (std::getline(in, line)) {
        std::stringstream ss { line };
        double val {};
        while (ss >> val) {
            array.push_back(val);
        }
        ++n_rows;
    }
    return { n_rows, static_cast<int>(array.size()) / n_rows };
}

// Read CKD for tests. Different CKD sections are read one by one from
// text files.
auto readCKD(const std::string& fixture_dir, tango::CKD& ckd) -> void
{
    // Work array for reading CKD contents
    std::vector<double> buf {};
    // Basic dimensions and the pixel mask
    const auto [n_rows, n_cols] { readArrayFromTxt(
      fixture_dir + "/ckd_pixel_mask.txt", buf) };
    ckd.n_detector_rows = n_rows;
    ckd.n_detector_cols = n_cols;
    ckd.npix = ckd.n_detector_rows * ckd.n_detector_cols;
    ckd.pixel_mask.resize(ckd.npix);
    for (int i {}; i < static_cast<int>(buf.size()); ++i) {
        ckd.pixel_mask[i] = static_cast<int>(std::round(buf[i])) == 1;
    }
    // Dark CKD
    ckd.dark.enabled = true;
    readArrayFromTxt(fixture_dir + "/ckd_dark_offset.txt", buf);
    ckd.dark.offset = buf;
    readArrayFromTxt(fixture_dir + "/ckd_dark_current.txt", buf);
    ckd.dark.current = buf;
    // Noise CKD
    ckd.noise.enabled = true;
    readArrayFromTxt(fixture_dir + "/ckd_conversion_gain.txt", buf);
    ckd.noise.g = buf;
    readArrayFromTxt(fixture_dir + "/ckd_read_noise.txt", buf);
    ckd.noise.n2 = buf;
    for (double& val : ckd.noise.n2) {
        val *= val;
    }
    // Nonlinearity CKD
    ckd.nonlin.enabled = true;
    readArrayFromTxt(fixture_dir + "/ckd_nonlin.txt", buf);
    std::vector<double> x_values(buf.size() / 2);
    std::vector<double> y_values(buf.size() / 2);
    for (int i {}; i < static_cast<int>(x_values.size()); ++i) {
        x_values[i] = buf[2 * i];
        y_values[i] = buf[2 * i + 1];
    }
    ckd.nonlin.spline = { x_values, y_values };
    // PRNU CKD
    ckd.prnu.enabled = true;
    readArrayFromTxt(fixture_dir + "/ckd_prnu.txt", buf);
    ckd.prnu.prnu = buf;
}

auto readL1(const std::string& fixture_dir, tango::L1& l1) -> void
{
    // Work array for reading L1 product contents
    std::vector<double> buf {};
    readArrayFromTxt(fixture_dir + "/detector_image.txt", buf);
    l1.image = buf;
    l1.pixel_mask.assign(buf.size(), false);
    l1.stdev.assign(buf.size(), 1.0);
    l1.exposure_time = 0.0460833;
    l1.nr_coadditions = 1;
}

auto absSum(std::vector<double> array) -> double
{
    return std::accumulate(
      array.begin(), array.end(), 0.0, [](const double s, const double a) {
          return s + std::abs(a);
      });
}
