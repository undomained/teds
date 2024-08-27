// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "ckd.h"

#include <cstdint>
#include <netcdf>
#include <spdlog/spdlog.h>

namespace tango {

// Retrieve values from a NetCDF variable into a buffer and reshape
// for destination.
static auto getAndReshape(const netCDF::NcVar& nc_var,
                          std::vector<std::vector<double>>& dest) -> void
{
    std::vector<double> buffer(dest.size() * dest.front().size());
    nc_var.getVar(buffer.data());
    for (int i {}; i < static_cast<int>(dest.size()); ++i) {
        for (int j {}; j < static_cast<int>(dest.at(i).size()); ++j) {
            dest[i][j] = buffer[i * dest.front().size() + j];
        }
    }
}

static auto getAndReshape3D(const netCDF::NcVar& nc_var,
                            std::vector<std::vector<std::vector<double>>>& dest)
  -> void
{
    int dim0 = dest.size();
    int dim1 = dest.front().size();
    int dim2 = dest.front().front().size();
    std::vector<double> buffer(dim0 * dim1 * dim2);
    nc_var.getVar(buffer.data());
    for (int i {}; i < dim0; ++i) {
        for (int j {}; j < dim1; ++j) {
            for (int k {}; k < dim2; ++k) {
                dest[i][j][k] = buffer[i * (dim1 * dim2) + j * dim2 + k];
            }
        }
    }
}

CKD::CKD(const std::string& filename, const double spectrum_width)
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    n_detector_rows = static_cast<int>(nc.getDim("detector_row").getSize());
    n_detector_cols = static_cast<int>(nc.getDim("detector_column").getSize());
    npix = n_detector_cols * n_detector_rows;
    n_act = static_cast<int>(nc.getDim("across_track").getSize());
    n_lbl = static_cast<int>(nc.getDim("lbl_samples").getSize());
    // Read the pixel mask which is possibly not yet in its final
    // state. It may be updated by one or more detector calibration
    // steps.
    std::vector<uint8_t> pixel_mask_u8(npix);
    nc.getVar("pixel_mask").getVar(pixel_mask_u8.data());
    pixel_mask.resize(npix);
    for (int i {}; i < npix; ++i) {
        pixel_mask[i] = static_cast<bool>(pixel_mask_u8[i]);
    }
    if (const netCDF::NcGroup grp { nc.getGroup("dark") }; !grp.isNull()) {
        spdlog::info("Reading dark CKD");
        dark.enabled = true;
        dark.offset.resize(npix);
        dark.current.resize(npix);
        grp.getVar("offset").getVar(dark.offset.data());
        grp.getVar("current").getVar(dark.current.data());
    }
    if (const netCDF::NcGroup grp { nc.getGroup("noise") }; !grp.isNull()) {
        spdlog::info("Reading noise CKD");
        noise.enabled = true;
        noise.g.resize(npix);
        noise.n2.resize(npix);
        grp.getVar("g").getVar(noise.g.data());
        grp.getVar("n").getVar(noise.n2.data());
        for (double& val : noise.n2) {
            val *= val;
        }
    }
    if (const netCDF::NcGroup grp { nc.getGroup("nonlinearity") };
        !grp.isNull()) {
        spdlog::info("Reading nonlinearity CKD");
        nonlin.enabled = true;
        const auto n_knots { grp.getDim("knots").getSize() };
        std::vector<double> knots(n_knots);
        std::vector<double> y(n_knots);
        grp.getVar("knots").getVar(knots.data());
        grp.getVar("y").getVar(y.data());
        nonlin.spline = { knots, y };
    }
    if (const netCDF::NcGroup grp { nc.getGroup("prnu") }; !grp.isNull()) {
        spdlog::info("Reading dark CKD");
        prnu.enabled = true;
        prnu.prnu.resize(npix);
        grp.getVar("prnu").getVar(prnu.prnu.data());
    }
    if (const netCDF::NcGroup grp { nc.getGroup("stray") }; !grp.isNull()) {
        spdlog::info("Reading stray light CKD");
        stray.enabled = true;
        stray.n_kernels = static_cast<int>(grp.getDim("kernel").getSize());
        stray.kernel_rows.resize(stray.n_kernels);
        grp.getVar("kernel_rows").getVar(stray.kernel_rows.data());
        stray.kernel_cols.resize(stray.n_kernels);
        grp.getVar("kernel_cols").getVar(stray.kernel_cols.data());
        stray.kernel_fft_sizes.resize(stray.n_kernels);
        grp.getVar("kernel_fft_sizes").getVar(stray.kernel_fft_sizes.data());
        stray.kernels_fft.resize(stray.n_kernels);
        for (int i {}; i < stray.n_kernels; ++i) {
            stray.kernels_fft[i].resize(stray.kernel_fft_sizes[i]);
        }
        std::vector<size_t> start { 0 };
        for (int i_kernel {}; i_kernel < stray.n_kernels; ++i_kernel) {
            // Need to multiply by 2 because of complex numbers
            const std::vector<size_t> count {
                2 * static_cast<size_t>(stray.kernel_fft_sizes.at(i_kernel))
            };
            grp.getVar("kernels_fft")
              .getVar(start, count, stray.kernels_fft.at(i_kernel).data());
            start.front() += count.front();
        }
        stray.eta.resize(npix);
        grp.getVar("eta").getVar(stray.eta.data());
        stray.weights.resize(stray.n_kernels, std::vector<double>(npix));
        std::vector<double> buf(stray.n_kernels * npix);
        grp.getVar("weights").getVar(buf.data());
        for (int i {}; i < stray.n_kernels; ++i) {
            for (int j {}; j < npix; ++j) {
                stray.weights[i][j] = buf[i * npix + j];
            }
        }
        stray.edges.resize(4 * stray.n_kernels);
        grp.getVar("edges").getVar(stray.edges.data());
    }
    if (const netCDF::NcGroup grp { nc.getGroup("swath") }; !grp.isNull()) {
        spdlog::info("Reading swath CKD");
        swath.enabled = true;
        swath.row_indices.resize(n_act, std::vector<double>(n_detector_cols));
        getAndReshape(grp.getVar("row_index"), swath.row_indices);
        genPixelIndices(spectrum_width);
    }
    if (const netCDF::NcGroup grp { nc.getGroup("spectral") }; !grp.isNull()) {
        spdlog::info("Reading spectral CKD");
        wave.enabled = true;
        wave.wavelength.resize(n_act, std::vector<double>(n_detector_cols));
        getAndReshape(grp.getVar("wavelength"), wave.wavelength);
    }
    if (const netCDF::NcGroup grp { nc.getGroup("radiometric") };
        !grp.isNull()) {
        spdlog::info("Reading radiometric CKD");
        rad.enabled = true;
        rad.rad.resize(n_detector_rows, std::vector<double>(n_detector_cols));
        getAndReshape(grp.getVar("radiometric"), rad.rad);
        // Resize 3D ISRF
        n_isrf_samples = static_cast<int>(grp.getDim("isrf_samples").getSize());
        rad.isrf.resize(n_act);
        rad.isrf_wl.resize(n_act);
        for (int i = 0; i < n_act; ++i) {
            rad.isrf[i].resize(n_lbl);
            rad.isrf_wl[i].resize(n_lbl);
            for (int j = 0; j < n_lbl; ++j) {
                rad.isrf[i][j].resize(n_isrf_samples);
                rad.isrf_wl[i][j].resize(n_isrf_samples);
            }
        }
        getAndReshape3D(grp.getVar("isrf"), rad.isrf);
        getAndReshape3D(grp.getVar("isrf_wavelengths"), rad.isrf_wl);
    }
}

auto CKD::genPixelIndices(const double spectrum_width) -> void
{
    // Convert fractional detector row coordinates into pixel indices
    // and weight factors.
    const double spectrum_width_half { spectrum_width / 2 };
    swath.n_indices = static_cast<int>(std::ceil(spectrum_width) + 1);
    swath.pix_indices.resize(
      n_act, std::vector<int>(n_detector_cols * swath.n_indices));
    swath.weights.resize(
      n_act, std::vector<double>(n_detector_cols * swath.n_indices));
    // Based on the argument spectrum_width, each spectral element
    // spans N rows on the detector. The weight of the first and last
    // rows are determined by where exactly the edge of the spectrum
    // falls. All middle rows have the same constant weight (1 before
    // normalization).
    for (int i_act {}; i_act < n_act; ++i_act) {
        for (int i {}; i < n_detector_cols; ++i) {
            const double row_first_d { swath.row_indices[i_act][i] + 0.5
                                       - spectrum_width_half };
            const double row_last_d { swath.row_indices[i_act][i] + 0.5
                                      + spectrum_width_half };
            const int row_first_i { static_cast<int>(row_first_d) };
            const int row_last_i { static_cast<int>(row_last_d) };
            auto& indices { swath.pix_indices[i_act] };
            auto& weights { swath.weights[i_act] };
            indices[i * swath.n_indices] = row_first_i * n_detector_cols + i;
            indices[(i + 1) * swath.n_indices - 1] =
              row_last_i * n_detector_cols + i;
            weights[i * swath.n_indices] = std::ceil(row_first_d) - row_first_d;
            weights[(i + 1) * swath.n_indices - 1] =
              row_last_d - std::floor(row_last_d);
            for (int i_ind { 1 }; i_ind < swath.n_indices - 1; ++i_ind) {
                const int idx { i * swath.n_indices + i_ind };
                indices[idx] = (row_first_i + i_ind) * n_detector_cols + i;
                weights[idx] = 1.0;
            }
            const int idx { i * swath.n_indices };
            if (indices[idx] == indices[idx + 1]) {
                weights[idx + 1] = 0.0;
            }
            if (indices[idx + swath.n_indices - 1]
                  == indices[idx + swath.n_indices - 2]
                && swath.n_indices > 2) {
                weights[idx + swath.n_indices - 2] = 0.0;
            }
            double weights_norm {};
            for (int i_ind {}; i_ind < swath.n_indices; ++i_ind) {
                weights_norm += weights[i * swath.n_indices + i_ind];
            }
            for (int i_ind {}; i_ind < swath.n_indices; ++i_ind) {
                weights[i * swath.n_indices + i_ind] /= weights_norm;
            }
        }
    }
}

} // namespace tango
