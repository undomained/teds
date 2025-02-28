#pragma once

#include <tango_l1b/ckd.h>
#include <tango_l1b/cubic_spline.h>
#include <tango_l1b/fourier.h>
#include <tango_l1b/l1.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <netcdf>
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

auto read2DArrayFromTxt(const std::string& filename,
                        std::vector<std::vector<double>>& array) -> void
{
    array.clear();
    std::ifstream in { filename };
    std::string line {};
    while (std::getline(in, line)) {
        std::stringstream ss { line };
        double val {};
        array.emplace_back(std::vector<double> {});
        while (ss >> val) {
            array.back().push_back(val);
        }
    }
}

// Because fixtures small return the flattened array directly
auto flatten2D(const std::vector<std::vector<double>>& data_2d)
  -> std::vector<double>
{
    const int n_rows { static_cast<int>(data_2d.size()) };
    const int n_cols { static_cast<int>(data_2d.front().size()) };
    std::vector<double> data_1d(n_rows * n_cols);
    for (int i {}; i < n_rows; ++i) {
        for (int j {}; j < n_cols; ++j) {
            data_1d[i * n_cols + j] = data_2d[i][j];
        }
    }
    return data_1d;
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
    // Detector mapping
    readArrayFromTxt(fixture_dir + "/ckd_swath_act_angle.txt",
                     ckd.swath.act_angles);
    readArrayFromTxt(fixture_dir + "/ckd_swath_wavelength.txt",
                     ckd.swath.wavelengths);
    readArrayFromTxt(fixture_dir + "/ckd_swath_act_map.txt", ckd.swath.act_map);
    readArrayFromTxt(fixture_dir + "/ckd_swath_wavelength_map.txt",
                     ckd.swath.wavelength_map);
    readArrayFromTxt(fixture_dir + "/ckd_swath_row_map.txt", ckd.swath.row_map);
    readArrayFromTxt(fixture_dir + "/ckd_swath_col_map.txt", ckd.swath.col_map);
    readArrayFromTxt(fixture_dir + "/ckd_swath_los.txt", ckd.swath.los);
    ckd.n_act = static_cast<int>(ckd.swath.act_angles.size());
    // Spectral
    read2DArrayFromTxt(fixture_dir + "/ckd_spectral.txt", ckd.wave.wavelengths);
    // Radiometric
    read2DArrayFromTxt(fixture_dir + "/ckd_radiometric.txt", ckd.rad.rad);
}

// Construct L1A product from fixtures
auto readL1A(const std::string& fixture_dir, tango::L1& l1) -> void
{
    l1.n_alt = 1;
    l1.exposure_time = 0.0460833;
    readArrayFromTxt(fixture_dir + "/detector_image.txt", l1.signal);
    l1.noise.assign(l1.signal.size(), 1.0);
    read2DArrayFromTxt(fixture_dir + "/l1b_wavelengths.txt", l1.wavelengths);
    readArrayFromTxt(fixture_dir + "/l1b_radiance.txt", l1.spectra);
    readArrayFromTxt(fixture_dir + "/l1b_radiance_stdev.txt", l1.spectra_noise);
}

// Read geolocation data from an L1B product
auto readGeo(const std::string& filename, tango::Geometry& geo) -> void
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    const auto n_alt { nc.getDim("along_track_sample").getSize() };
    const auto n_act { nc.getDim("across_track_sample").getSize() };
    const auto n_bins { n_alt * n_act };
    geo.lat.resize(n_bins);
    geo.lon.resize(n_bins);
    geo.height.resize(n_bins);
    geo.sza.resize(n_bins);
    geo.saa.resize(n_bins);
    geo.vza.resize(n_bins);
    geo.vaa.resize(n_bins);
    const auto grp { nc.getGroup("geolocation_data") };
    grp.getVar("latitude").getVar(geo.lat.data());
    grp.getVar("longitude").getVar(geo.lon.data());
    grp.getVar("height").getVar(geo.height.data());
    grp.getVar("sensor_zenith").getVar(geo.vza.data());
    grp.getVar("sensor_azimuth").getVar(geo.vaa.data());
    grp.getVar("solar_zenith").getVar(geo.sza.data());
    grp.getVar("solar_azimuth").getVar(geo.saa.data());
}

// Write CKD contents to a NetCDF file so that the read function can
// be properly tested.
auto writeCKD(const std::string& fixture_dir,
              const std::string& filename,
              const tango::CKD& ckd) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    // Dimensions
    const auto nc_n_rows { nc.addDim("detector_row", ckd.n_detector_rows) };
    const auto nc_n_cols { nc.addDim("detector_column", ckd.n_detector_cols) };
    const auto nc_n_act { nc.addDim("across_track_sample", ckd.n_act) };
    const auto nc_vector { nc.addDim("vector", 3) };
    const std::vector<netCDF::NcDim> nc_detector_shape { nc_n_rows, nc_n_cols };
    netCDF::NcGroup nc_grp {};
    netCDF::NcVar nc_var {};
    // Pixel mask
    std::vector<uint8_t> pixel_mask_u8(ckd.npix);
    for (int i {}; i < ckd.npix; ++i) {
        pixel_mask_u8[i] = static_cast<uint8_t>(ckd.pixel_mask[i]);
    }
    nc.addVar("pixel_mask", netCDF::ncByte, nc_detector_shape)
      .putVar(pixel_mask_u8.data());
    // Dark CKD
    nc_grp = nc.addGroup("dark");
    nc_grp.addVar("offset", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.dark.offset.data());
    nc_grp.addVar("current", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.dark.current.data());
    // Noise CKD
    nc_grp = nc.addGroup("noise");
    nc_grp.addVar("g", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.noise.g.data());
    std::vector<double> buf { ckd.noise.n2 };
    for (double& val : buf) {
        val = std::sqrt(val);
    }
    nc_grp.addVar("n", netCDF::ncDouble, nc_detector_shape).putVar(buf.data());
    // Nonlinearity CKD
    std::vector<double> x_values {};
    std::vector<double> y_values {};
    readNonlinCKD(fixture_dir, x_values, y_values);
    nc_grp = nc.addGroup("nonlinearity");
    const auto nc_knots { nc_grp.addDim("knots", x_values.size()) };
    nc_grp.addVar("knots", netCDF::ncDouble, nc_knots).putVar(x_values.data());
    nc_grp.addVar("y", netCDF::ncDouble, nc_knots).putVar(y_values.data());
    // PRNU CKD
    nc.addGroup("prnu")
      .addVar("prnu", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.prnu.prnu.data());
    // Stray light CKD
    nc_grp = nc.addGroup("stray");
    const auto nc_n_kernel { nc_grp.addDim("kernel", ckd.stray.n_kernels) };
    const auto& _s { ckd.stray.kernel_fft_sizes };
    const auto kernel_fft_size { 2
                                 * std::accumulate(_s.cbegin(), _s.cend(), 0) };
    const auto nc_kernels_fft_size { nc_grp.addDim("kernels_fft_size",
                                                   kernel_fft_size) };
    const auto nc_edges_of_box { nc_grp.addDim("edges_of_box", tango::box::n) };
    nc_grp.addVar("kernel_rows", netCDF::ncInt, nc_n_kernel)
      .putVar(ckd.stray.kernel_rows.data());
    nc_grp.addVar("kernel_cols", netCDF::ncInt, nc_n_kernel)
      .putVar(ckd.stray.kernel_cols.data());
    nc_grp.addVar("kernel_fft_sizes", netCDF::ncInt, nc_n_kernel)
      .putVar(ckd.stray.kernel_fft_sizes.data());
    auto var { nc_grp.addVar(
      "kernels_fft", netCDF::ncDouble, nc_kernels_fft_size) };
    std::vector<size_t> start { 0 };
    for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
        std::vector<size_t> count {
            2 * static_cast<size_t>(ckd.stray.kernel_fft_sizes[i_kernel])
        };
        var.putVar(start, count, ckd.stray.kernels_fft[i_kernel].data());
        start.front() += count.front();
    }
    nc_grp.addVar("eta", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.stray.eta.data());
    nc_var = nc_grp.addVar(
      "weights", netCDF::ncDouble, { nc_n_kernel, nc_n_rows, nc_n_cols });
    buf.resize(ckd.stray.n_kernels * ckd.npix);
    for (int i {}; i < ckd.stray.n_kernels; ++i) {
        for (int j {}; j < ckd.npix; ++j) {
            buf[i * ckd.npix + j] = ckd.stray.weights[i][j];
        }
    }
    nc_var.putVar(buf.data());
    nc_grp.addVar("edges", netCDF::ncInt, { nc_n_kernel, nc_edges_of_box })
      .putVar(ckd.stray.edges.data());
    // Swath CKD
    nc_grp = nc.addGroup("swath");
    const auto nc_swath_n_waves { nc_grp.addDim("wavelength",
                                                ckd.swath.wavelengths.size()) };
    nc_grp.addVar("wavelength", netCDF::ncDouble, nc_swath_n_waves)
      .putVar(ckd.swath.wavelengths.data());
    nc_grp.addVar("act_angle", netCDF::ncDouble, nc_n_act)
      .putVar(ckd.swath.act_angles.data());
    nc_grp.addVar("act_map", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.swath.act_map.data());
    nc_grp.addVar("wavelength_map", netCDF::ncDouble, nc_detector_shape)
      .putVar(ckd.swath.wavelength_map.data());
    nc_grp.addVar("row_map", netCDF::ncDouble, { nc_n_act, nc_swath_n_waves })
      .putVar(ckd.swath.row_map.data());
    nc_grp.addVar("col_map", netCDF::ncDouble, { nc_n_act, nc_swath_n_waves })
      .putVar(ckd.swath.col_map.data());
    nc_grp.addVar("line_of_sight", netCDF::ncDouble, { nc_n_act, nc_vector })
      .putVar(ckd.swath.los.data());
    // Spectral CKD
    const auto wavelengths { flatten2D(ckd.wave.wavelengths) };
    nc.addGroup("spectral")
      .addVar("wavelength", netCDF::ncDouble, { nc_n_act, nc_n_cols })
      .putVar(wavelengths.data());
    // Radiometric CKD
    const auto rad { flatten2D(ckd.rad.rad) };
    nc.addGroup("radiometric")
      .addVar("radiometric", netCDF::ncDouble, { nc_n_act, nc_n_cols })
      .putVar(rad.data());
}

// Write L1A product to NetCDF file with minimal structure (no
// unnecessary attributes).
auto writeL1A(const std::string& fixture_dir,
              const std::string& filename,
              tango::L1& l1) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    // Dimensions
    const auto nc_n_alt { nc.addDim("along_track_sample", l1.n_alt) };
    const auto nc_n_bin { nc.addDim("bin", l1.signal.size() / l1.n_alt) };
    nc.putAtt("product_type", "L1A");
    nc.putAtt("instrument", "Carbon");
    // Image attributes
    auto grp { nc.addGroup("image_attributes") };
    l1.time = { 44025.4149999619 };
    grp.addVar("time", netCDF::ncDouble, nc_n_alt).putVar(l1.time.data());
    l1.tai_seconds = { 2038997625 };
    grp.addVar("tai_seconds", netCDF::ncUint, nc_n_alt)
      .putVar(l1.tai_seconds.data());
    l1.tai_subsec = { 27197 };
    grp.addVar("tai_subsec", netCDF::ncUint, nc_n_alt)
      .putVar(l1.tai_subsec.data());
    grp.addVar("binning_table", netCDF::ncByte).putVar(&l1.binning_table_id);
    l1.nr_coadditions = 2;
    grp.addVar("nr_coadditions", netCDF::ncUshort).putVar(&l1.nr_coadditions);
    l1.exposure_time = 0.01724385;
    grp.addVar("exposure_time", netCDF::ncDouble).putVar(&l1.exposure_time);
    // Science data
    nc.addGroup("science_data")
      .addVar("detector_image", netCDF::ncInt, { nc_n_alt, nc_n_bin })
      .putVar(l1.signal.data());
    // Navigation data
    std::vector<double> nav_times {};
    readArrayFromTxt(fixture_dir + "/navigation_time.txt", nav_times);
    readArrayFromTxt(fixture_dir + "/navigation_orb_pos.txt", l1.orb_pos);
    std::vector<double> att_quat {};
    readArrayFromTxt(fixture_dir + "/navigation_att_quat.txt", att_quat);
    grp = nc.addGroup("navigation_data");
    const auto dim_time { grp.addDim("time", nav_times.size()) };
    const auto dim_vec { grp.addDim("vector_elements", 3) };
    const auto dim_quat { grp.addDim("quaternion_elements", 4) };
    grp.addVar("time", netCDF::ncDouble, dim_time).putVar(nav_times.data());
    grp.addVar("orb_pos", netCDF::ncDouble, { dim_time, dim_vec })
      .putVar(l1.orb_pos.data());
    grp.addVar("att_quat", netCDF::ncDouble, { dim_time, dim_quat })
      .putVar(att_quat.data());
}

// Write SGM product to NetCDF file
auto writeSGM(const std::string& fixture_dir,
              const std::string& filename,
              const tango::CKD& ckd,
              tango::L1& l1) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    // Line-by-line (LBL) wavelength grid is based on the CKD
    // intermediate wavelength grid.
    constexpr int lbl_wave_multiplier { 3 };
    const auto nc_n_waves { nc.addDim(
      "wavelength", lbl_wave_multiplier * ckd.swath.wavelengths.size()) };
    const auto nc_n_alt { nc.addDim("along_track_sample", l1.n_alt) };
    const auto nc_n_act { nc.addDim("across_track_sample", ckd.n_act) };
    nc.putAtt("product_type", "SGM");
    // Generate LBL grid
    std::vector<double> lbl_wavelengths(nc_n_waves.getSize());
    for (int i {}; i < static_cast<int>(lbl_wavelengths.size()); ++i) {
        const double& wave_0 { ckd.swath.wavelengths.front() };
        const double& wave_n { ckd.swath.wavelengths.back() };
        constexpr double margin { 1.0 };
        lbl_wavelengths[i] =
          wave_0 - margin
          + i * (wave_n - wave_0 + 2 * margin) / (lbl_wavelengths.size() - 1);
    }
    // Generate LBL spectra
    std::vector<double> lbl_spectra(ckd.n_act * lbl_wavelengths.size());
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        const tango::CubicSpline spline {
            ckd.wave.wavelengths[i_act],
            { l1.spectra.begin() + i_act * ckd.n_detector_cols,
              l1.spectra.begin() + (i_act + 1) * ckd.n_detector_cols }
        };
        for (int i {}; i < static_cast<int>(lbl_wavelengths.size()); ++i) {
            lbl_spectra[i_act * lbl_wavelengths.size() + i] =
              spline.eval(lbl_wavelengths[i]);
        }
    }
    // Write to NetCDF file
    nc.addVar("wavelength", netCDF::ncDouble, nc_n_waves)
      .putVar(lbl_wavelengths.data());
    nc.addVar("radiance", netCDF::ncDouble, { nc_n_alt, nc_n_act, nc_n_waves })
      .putVar(lbl_spectra.data());
}

auto writeGeometry(const std::string& filename, tango::L1& l1) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    // Dimensions
    const auto nc_n_alt { nc.addDim("along_track_sample", l1.n_alt) };
    l1.time = { 44025.4149999619 };
    auto nc_time { nc.addVar("time", netCDF::ncDouble, nc_n_alt) };
    nc_time.putVar(l1.time.data());
    nc_time.putAtt("units", "seconds since 2022-08-12");
    l1.tai_seconds = { 2038997625 };
    nc.addVar("tai_seconds", netCDF::ncUint, nc_n_alt)
      .putVar(l1.tai_seconds.data());
    l1.tai_subsec = { 27197 };
    nc.addVar("tai_subsec", netCDF::ncUint, nc_n_alt)
      .putVar(l1.tai_subsec.data());
}

auto writeNavigation(const std::string& fixture_dir,
                     const std::string& filename) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    std::vector<double> buf {};
    readArrayFromTxt(fixture_dir + "/navigation_time.txt", buf);
    const auto dim_time { nc.addDim("time", buf.size()) };
    const auto dim_vec { nc.addDim("vector_elements", 3) };
    const auto dim_quat { nc.addDim("quaternion_elements", 4) };
    nc.addVar("time", netCDF::ncDouble, dim_time).putVar(buf.data());
    readArrayFromTxt(fixture_dir + "/navigation_orb_pos.txt", buf);
    nc.addVar("orb_pos", netCDF::ncDouble, { dim_time, dim_vec })
      .putVar(buf.data());
    readArrayFromTxt(fixture_dir + "/navigation_att_quat.txt", buf);
    nc.addVar("att_quat", netCDF::ncDouble, { dim_time, dim_quat })
      .putVar(buf.data());
}

// Create binning tables and write to NetCDF file
auto writeBinningTable(const std::string& fixture_dir,
                       const std::string& filename) -> void
{
    std::vector<double> buf {};
    const auto [n_rows, n_cols] { readArrayFromTxt(
      fixture_dir + "/detector_image.txt", buf) };
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    const auto nc_rows { nc.addDim("row", n_rows) };
    const auto nc_cols { nc.addDim("column", n_cols) };
    for (int bin { 1 }; bin <= 5; ++bin) {
        // Binning table
        std::vector<uint32_t> binning_table(n_rows * n_cols, 0);
        const int n_rows_red { n_rows / bin };
        int n_bins { n_rows_red * n_cols };
        uint32_t idx {};
        for (int i {}; i < n_rows_red; ++i) {
            for (int j {}; j < n_cols; ++j) {
                for (int k {}; k < bin; ++k) {
                    binning_table[(i * bin + k) * n_cols + j] = idx;
                }
                ++idx;
            }
        }
        const int n_rows_remaining { n_rows - n_rows_red * bin };
        if (n_rows_remaining > 0) {
            n_bins += n_cols;
        }
        const int bin_new { std::max(1, n_rows_remaining) };
        // Append to binning table
        if (n_rows_remaining > 0) {
            for (int j {}; j < n_cols; ++j) {
                for (int k {}; k < bin_new; ++k) {
                    binning_table[(n_rows - n_rows_remaining + k) * n_cols
                                  + j] = idx;
                }
                ++idx;
            }
        }
        // Count table
        std::vector<uint16_t> count_table(n_bins, bin);
        if (n_rows_remaining > 0) {
            for (int j {}; j < n_cols; ++j) {
                count_table[count_table.size() - n_cols + j] = bin_new;
            }
        }
        // Write to NetCDF file
        auto grp { nc.addGroup("Table_" + std::to_string(bin)) };
        const auto nc_bins { grp.addDim("bins", n_bins) };
        grp.addVar("binning_table", netCDF::ncUint, { nc_rows, nc_cols })
          .putVar(binning_table.data());
        grp.addVar("count_table", netCDF::ncUshort, nc_bins)
          .putVar(count_table.data());
    }
}

auto absSum(std::vector<double> array) -> double
{
    return std::accumulate(
      array.begin(), array.end(), 0.0, [](const double s, const double a) {
          return s + std::abs(a);
      });
}
