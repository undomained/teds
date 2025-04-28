// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "ckd.h"

#include "binning_table.h"

#include <netcdf>
#include <spdlog/spdlog.h>

namespace tango {

CKD::CKD(const std::string& filename)
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    n_detector_rows = static_cast<int>(nc.getDim("detector_row").getSize());
    n_detector_cols = static_cast<int>(nc.getDim("detector_column").getSize());
    npix = n_detector_cols * n_detector_rows;
    n_detector_rows_binned = n_detector_rows;
    n_detector_cols_binned = n_detector_cols;
    npix_binned = npix;
    n_act = static_cast<int>(nc.getDim("across_track_sample").getSize());
    const auto n_wavelengths { nc.getDim("wavelength").getSize() };
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
    Eigen::ArrayXd knots(n_knots);
    Eigen::ArrayXd y(n_knots);
    grp.getVar("observed").getVar(knots.data());
    grp.getVar("expected").getVar(y.data());
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
        stray.kernels_fft[i].resize(stray.kernel_fft_sizes(i));
    }
    std::vector<size_t> start { 0 };
    for (int i_kernel {}; i_kernel < stray.n_kernels; ++i_kernel) {
        // Need to multiply by 2 because of complex numbers
        const std::vector<size_t> count {
            2 * static_cast<size_t>(stray.kernel_fft_sizes(i_kernel))
        };
        grp.getVar("kernels_fft")
          .getVar(start, count, stray.kernels_fft.at(i_kernel).data());
        start.front() += count.front();
    }
    stray.eta.resize(n_detector_rows, n_detector_cols);
    grp.getVar("eta").getVar(stray.eta.data());
    stray.weights.resize(stray.n_kernels, npix);
    grp.getVar("weights").getVar(stray.weights.data());
    stray.edges.resize(stray.n_kernels, box::n);
    grp.getVar("edges").getVar(stray.edges.data());

    spdlog::info("Reading swath CKD");
    grp = nc.getGroup("swath");
    swath.act_angles.resize(n_act);
    grp.getVar("act_angle").getVar(swath.act_angles.data());
    swath.act_map.resize(n_detector_rows, n_detector_cols);
    grp.getVar("act_map").getVar(swath.act_map.data());
    swath.wavelength_map.resize(n_detector_rows, n_detector_cols);
    grp.getVar("wavelength_map").getVar(swath.wavelength_map.data());
    swath.row_map.resize(n_act, n_wavelengths);
    grp.getVar("row_map").getVar(swath.row_map.data());
    swath.col_map.resize(n_act, n_wavelengths);
    grp.getVar("col_map").getVar(swath.col_map.data());
    swath.los.resize(n_act, dims::vec);
    grp.getVar("line_of_sight").getVar(swath.los.data());

    spdlog::info("Reading spectral CKD");
    grp = nc.getGroup("spectral");
    wave.wavelengths.resize(n_wavelengths);
    grp.getVar("wavelength").getVar(wave.wavelengths.data());

    spdlog::info("Reading radiometric CKD");
    grp = nc.getGroup("radiometric");
    rad.rad.resize(n_act, n_wavelengths);
    grp.getVar("radiometric").getVar(rad.rad.data());

    // Set CKD values of bad detector pixels to 1 for numerical stability.
    for (int i {}; i < npix; ++i) {
        if (pixel_mask(i)) {
            dark.offset(i) = 1.0;
            dark.current(i) = 1.0;
            noise.g(i) = 1.0;
            noise.n2(i) = 1.0;
            prnu.prnu(i) = 1.0;
        }
    }
}

auto CKD::bin(const BinningTable& binning_table) -> void
{
    pixel_mask = binning_table.bin(pixel_mask);
    dark.offset = binning_table.bin(dark.offset);
    dark.current = binning_table.bin(dark.current);
    noise.g = binning_table.bin(noise.g);
    noise.n2 = binning_table.bin(noise.n2);
    prnu.prnu = binning_table.bin(prnu.prnu);
    // For now, assume no binning across columns and thus binned and
    // unbinned number of columns are always the same.
    npix_binned = binning_table.nBins();
    n_detector_rows_binned = npix_binned / n_detector_cols_binned;
}

} // namespace tango
