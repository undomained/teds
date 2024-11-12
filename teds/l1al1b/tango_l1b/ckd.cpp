// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "ckd.h"

#include "constants.h"

#include <algorithm>
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

CKD::CKD(const std::string& filename)
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    n_detector_rows = static_cast<int>(nc.getDim("detector_row").getSize());
    n_detector_cols = static_cast<int>(nc.getDim("detector_column").getSize());
    npix = n_detector_cols * n_detector_rows;
    n_act = static_cast<int>(nc.getDim("across_track_sample").getSize());
    // Read the pixel mask which is possibly not yet in its final
    // state. It may be updated by one or more detector calibration
    // steps.
    std::vector<uint8_t> pixel_mask_u8(npix);
    nc.getVar("pixel_mask").getVar(pixel_mask_u8.data());
    pixel_mask.resize(npix);
    for (int i {}; i < npix; ++i) {
        pixel_mask[i] = static_cast<bool>(pixel_mask_u8[i]);
    }

    spdlog::info("Reading dark CKD");
    netCDF::NcGroup grp {};
    grp = nc.getGroup("dark");
    dark.offset.resize(npix);
    dark.current.resize(npix);
    grp.getVar("offset").getVar(dark.offset.data());
    grp.getVar("current").getVar(dark.current.data());

    spdlog::info("Reading noise CKD");
    grp = nc.getGroup("noise");
    noise.g.resize(npix);
    noise.n2.resize(npix);
    grp.getVar("g").getVar(noise.g.data());
    grp.getVar("n").getVar(noise.n2.data());
    for (double& val : noise.n2) {
        val *= val;
    }

    spdlog::info("Reading nonlinearity CKD");
    grp = nc.getGroup("nonlinearity");
    const auto n_knots { grp.getDim("knots").getSize() };
    std::vector<double> knots(n_knots);
    std::vector<double> y(n_knots);
    grp.getVar("knots").getVar(knots.data());
    grp.getVar("y").getVar(y.data());
    nonlin.spline = { knots, y };

    spdlog::info("Reading dark CKD");
    grp = nc.getGroup("prnu");
    prnu.prnu.resize(npix);
    grp.getVar("prnu").getVar(prnu.prnu.data());

    spdlog::info("Reading stray light CKD");
    grp = nc.getGroup("stray");
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
    std::vector<double> buf(stray.n_kernels* npix);
    grp.getVar("weights").getVar(buf.data());
    for (int i {}; i < stray.n_kernels; ++i) {
        for (int j {}; j < npix; ++j) {
            stray.weights[i][j] = buf[i * npix + j];
        }
    }
    stray.edges.resize(4 * stray.n_kernels);
    grp.getVar("edges").getVar(stray.edges.data());

    spdlog::info("Reading swath CKD");
    grp = nc.getGroup("swath");
    const auto n_intermediate_waves { grp.getDim("wavelength").getSize() };
    swath.wavelengths.resize(n_intermediate_waves);
    grp.getVar("wavelength").getVar(swath.wavelengths.data());
    swath.act_angles.resize(n_act);
    grp.getVar("act_angle").getVar(swath.act_angles.data());
    swath.act_map.resize(npix);
    grp.getVar("act_map").getVar(swath.act_map.data());
    swath.wavelength_map.resize(npix);
    grp.getVar("wavelength_map").getVar(swath.wavelength_map.data());
    swath.row_map.resize(n_act * n_intermediate_waves);
    grp.getVar("row_map").getVar(swath.row_map.data());
    swath.col_map.resize(n_act * n_intermediate_waves);
    grp.getVar("col_map").getVar(swath.col_map.data());
    swath.los.resize(n_act * dims::vec);
    grp.getVar("line_of_sight").getVar(swath.los.data());

    spdlog::info("Reading spectral CKD");
    grp = nc.getGroup("spectral");
    wave.wavelengths.resize(n_act, std::vector<double>(n_detector_cols));
    getAndReshape(grp.getVar("wavelength"), wave.wavelengths);
    if (wave.wavelengths.front()[1] < wave.wavelengths.front()[0]) {
        for (auto& waves : wave.wavelengths) {
            std::ranges::reverse(waves);
        }
    }

    spdlog::info("Reading radiometric CKD");
    grp = nc.getGroup("radiometric");
    rad.rad.resize(n_act, std::vector<double>(n_detector_cols));
    getAndReshape(grp.getVar("radiometric"), rad.rad);
}

} // namespace tango
