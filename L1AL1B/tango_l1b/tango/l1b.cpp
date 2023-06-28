// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "array.h"
#include "matrix.h"
#include "lininv.h"
#include "vector.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "planet.h"
#include "ckd.h"
#include "binningtable.h"
#include "l1a.h"
#include "l1b.h"

#include <algorithm>
#include <cstring>

double AngToX(const double ang, const double h, const double eR){
	double xl;
	//const double eR = 6371.0;
	const long double pi = 3.141592653589793;
	const double a = 1;
	const double b = -2.0 * (eR + h) * cos(ang * 2 * pi / 360);
	const double c = pow((eR + h), 2) - pow(eR, 2);
	const double l = (-b - pow((pow(b,2)-4 * a * c),0.5)) / (2.0 * a);
	const double fi = acos((pow(eR, 2) + pow((eR + h), 2) - pow(l, 2)) / (2 * eR * (eR + h)));

	xl = eR * fi;
	if(ang < 0){
		xl = -xl;
	}

	return xl;
}

// Settings functions.
Settings_l1b::Settings_l1b( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "l1b";
    l1a_type = L1A_FLIGHT; // Default L1A style.
    overall_percentages = true; // Turn on overall progress bar.
} // }}}

int Settings_l1b::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(outputfile); // Output L1B file.
    recognize_setting(order); // B-spline order of wavelength dependence of q and u.
    recognize_setting_vector(knots); // Wavelengths of knot positions in B-spline.
    recognize_setting(wave_min); // Wavelength to start the demodulation window.
    recognize_setting(wave_max); // Wavelength to end the demodulation window.
    recognize_setting_vector(l1b_wavelength); // Central wavelengths of polarimetric product.
    recognize_setting_vector(intensity_wavelength); // Wavelengths for intensity-only grid.
    recognize_setting(resolving_power); // Full-width half-maximum of Gauss convolution of polarimetric L1B product.
    recognize_setting(gauss_range); // Range of Gauss in FWHMs on both sides.
    recognize_setting(demod_noise); // Flag to also calculate noise propagation on demodulation products.
    recognize_setting(geolocation); // Execution flag for geolocation (Turn off when using OCAL measurements as test).
    recognize_setting(utcfile); // UTC Time-difference file.
    recognize_setting(demfile); // Detailed elevation map. Leave empty for placeholder DEM function.
    recognize_setting(semi_major_axis); // Semi-major axis of the earth.
    recognize_setting(semi_minor_axis); // Semi-minor axis of the earth.
    recognize_setting(latitude_tol); // Tolerance for iterative calculation of latitude during geolocation (radians).
    recognize_setting(mountain); // Highest mountain to expect.
    recognize_setting(movedistance); // Distance to move in one iteration avoiding to skip mountains.
    recognize_setting(extremeweight); // Clip on weight factor for guess during geolocation convergence.
    recognize_setting(geolocation_tol); // Tolerance for convergence along line of sight (length units).
    recognize_setting(l1b_output_level); // Level of detailed output. Own implementation, because output file is not the CKD file.
    recognize_setting(first_proc_rel_workload);
    recognize_setting(mpi_send_size);
    recognize_setting(geometry_file);

    return 0;

} // }}}

// Constructor for the L1B CKD structure.
L1B::L1B( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("l1b");
    set = make_unique<Settings_l1b>(this);
    Processor::set = set.get();
} // }}}
L1B::~L1B() {}

// Initialization before pixel loop.
int L1B::process_init( // {{{
)
{

    // Some administration to begin with.
    // If you skip the swath part of your CKD, the geolocation will be
    // skipped. Because there may be other reasons to skip the geolocation,
    // we pool to the one custom flag.
    if (ckd->lev <= LEVEL_SWATHCAL) set->geolocation = false;
    if (set->geolocation && ckd->swath_skip) {
        writelog(log_warning,"Warning: Swath vector step was skipped. No geolocation can be calculated.");
        set->geolocation = false;
    }

    // Check legal settings.
    check_error(set->outputfile.compare("") == 0,"Error: Missing output file name. Setting 'outputfile'");
    if (set->geolocation) {
        /*
        check_error(set->utcfile.compare("") == 0,"Error: Missing input UTC file. Setting 'utcfile'");
        check_error(set->semi_major_axis == NC_FILL_DOUBLE,"Error: Missing semi-major axis of the planet. Setting 'semi_major_axis'");
        check_error(set->semi_minor_axis == NC_FILL_DOUBLE,"Error: Missing semi-minor axis of the planet. Setting 'semi_minor_axis'");
        check_error(set->latitude_tol == NC_FILL_DOUBLE,"Error: Missing tolerance for iterative calculation of latitude during geolocation (radians). Setting 'latitude_tol'");
        check_error(set->mountain == NC_FILL_DOUBLE,"Error: Missing estimate of highest mountain (any higher point is always considered to be above the surface). Setting 'mountain'");
        check_error(set->movedistance == NC_FILL_DOUBLE,"Error: Missing sampling distance for scanning line of sight.");
        check_error(set->extremeweight == NC_FILL_DOUBLE,"Error: Missing extreme-weight parameter for optimizing geo-location. Any value between zero and a half make sense (not zero itself). Setting 'extremeweight'");
        check_error(set->geolocation_tol == NC_FILL_DOUBLE,"Error: Missing tolerance for convergence along line of sight (length units). Setting 'geolocation_tol'");
        */
    }

    // Because of repeated if-clauses, declarations are done before.
    // Define dimensions.
    NcDim dimid_obs; // Number of images.
    NcDim dimid_fov; // Number of swath positions per viewport.
    NcDim dimid_int_wave; // Number of intensity (raw) wavelengths.
    NcDim dimid_pol_wave; // Number of L1B output wavelengths.

    nc_l1b = make_unique<NetCDF_object>(this);
    handle(nc_l1b->open(set->outputfile, NcFile::replace));

    // Determine how many images each MPI process needs to handle
    std::vector<int> workloads(num_procs);
    // In a serial run the relative load on process 0 cannot be variable.
    if (num_procs == 1) {
        set->first_proc_rel_workload = 1.0;
    }
    // The first one is assigned a fixed number determined by an input setting.
    workloads.front() = static_cast<int>(std::round(set->first_proc_rel_workload
                                                    * nl1a_total / num_procs));
    // Distribute the remaining images equally among the remaining MPI
    // processes.
    const int nl1a_remaining { static_cast<int>(nl1a_total)
                               - workloads.front() };
    for (int i_rank { 1 }; i_rank < static_cast<int>(workloads.size());
         ++i_rank) {
        workloads.at(i_rank) =
          static_cast<int>(nl1a_remaining) / (num_procs - 1);
    }
    // Some processes might be assigned one more image than others.
    if (num_procs > 1) {
        for (int i { 1 };
             i < static_cast<int>(nl1a_remaining) % (num_procs - 1) + 1;
             ++i) {
            ++workloads.at(i);
        }
    }
    mpi_ranges.resize(num_procs + 1);
    mpi_ranges.front() = 0;
    for (int i_rank { 1 }; i_rank < num_procs + 1; ++i_rank) {
        mpi_ranges.at(i_rank) =
          mpi_ranges.at(i_rank - 1) + workloads.at(i_rank - 1);
    }
    if (num_procs > 1) {
        if (num_procs > nl1a_total) {
            writelog(log_error,
                     "Error: number of MPI tasks should not be less than the "
                     "number of images.");
            return 1;
        }
        writelog(log_info,
                 "Images are distributed among the MPI processes as follows:");
        std::stringstream workloads_line {};
        for (const auto& el : workloads) {
            workloads_line << "  " << el;
        }
        writelog(log_info, "%s", workloads_line.str().c_str());
    }

    // Define L1B-specific global attributes.
    // TODO: Add CDL version date when CDL is finished.
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("cdm_data_type","swath"));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("processing_lsciel","L1B"));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("title","TANGO Carbon Level-1B Data"));
    // Write global attributes from orbit data.
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("startDirection",orb.startDirection));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("endDirection",orb.endDirection));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("orbit_number",ncInt64,orb.orbit_number));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("time_coverage_start",orb.time_coverage_start));
    netcdf_check(nc_l1b,nc_l1b->ncid->putAtt("time_coverage_stop",orb.time_coverage_stop));

    // Define groups.
    NcGroup grp_bin;
    NcGroup grp_geolocation;
    NcGroup grp_observation;
    netcdf_check(nc_l1b,grp_bin = nc_l1b->ncid->addGroup("IMAGE_ATTRIBUTES"));
    netcdf_check(nc_l1b,grp_geolocation = nc_l1b->ncid->addGroup("GEOLOCATION_DATA"));
    netcdf_check(nc_l1b,grp_observation = nc_l1b->ncid->addGroup("OBSERVATION_DATA"));
    NcVar var; // Temporary variable identifier for variables that are directly filled.

    // Define dimensions.
    netcdf_check(nc_l1b,dimid_obs = nc_l1b->ncid->addDim("bins_along_track",nl1a_total));

    // Define variables.
    // Variables will be defined if the calculation reaches that point. Operationally, all
    // calculations will be done, so all variables will be created. For testing, it is possible
    // to just calculate radiance, or not even that. Also, geolocation can be skipped.

    // The image time is just a copy from the L1A, so always exists.
    netcdf_check(nc_l1b,var_time = grp_bin.addVar("image_time",ncDouble,dimid_obs));
    netcdf_check(nc_l1b,var_time.putAtt("long_name","image time (seconds of day)"));
    netcdf_check(nc_l1b,var_time.putAtt("valid_min",ncDouble,0.0));
    netcdf_check(nc_l1b,var_time.putAtt("valid_max",ncDouble,86400.999999));
    netcdf_check(nc_l1b,var_time.putAtt("units","seconds"));
    netcdf_check(nc_l1b,var_time.putAtt("reference",orb.time_reference));

    // If anything is calculated related to a spectrum. This can be radiance, demodulation or geolocation.
    if (ckd->lev > LEVEL_RADCAL || set->geolocation) {
        netcdf_check(nc_l1b,dimid_fov = nc_l1b->ncid->addDim("bins_across_track",ckd->dim_fov));
        // Map of viewport indices per FOV.
        vector<uint8_t> ivp_fov(ckd->dim_fov);
        size_t ivp = 0;
        for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
            ivp_fov[ifov] = (uint8_t) ivp;
        }
    }

    // Everything for radiance.
    if (ckd->lev > LEVEL_RADCAL) {

        // Define L1B wavelengths for intensity grid. This will be updated
        // to the raw grid with binning. TODO: Do this.
        // Now, we read intensity wavelengths from the settings.
        dim_int_wave = ckd->dim_detector_spec;

        // The dimension for intensity wavelengths.
        netcdf_check(nc_l1b,dimid_int_wave = nc_l1b->ncid->addDim("bins_spectral",dim_int_wave));

        vector<NcDim> dims_int = {dimid_int_wave};

        // Wavelengths intensity grid.
        var = grp_observation.addVar(
          "wavelengths", ncDouble, { dimid_fov, dimid_int_wave });
        var.putAtt("long_name", "wavelength at intensity samples");
        var.putAtt("units", "nm");
        vector<double> wavelengths_reversed { ckd->wave_target };
        if (ckd->wave_target.front() > ckd->wave_target.back()) {
            for (int i {}; i < ckd->dim_fov; ++i) {
                for (int j {}; j < ckd->dim_detector_spec / 2; ++j) {
                    std::swap(
                      wavelengths_reversed[i * ckd->dim_detector_spec + j],
                      wavelengths_reversed[i * ckd->dim_detector_spec
                                           + ckd->dim_detector_spec - 1 - j]);
                }
            }
        }
        var.putVar(wavelengths_reversed.data());

        // The wavelengths and the intensity-only spectra.
        vector<NcDim> dims_spectra_int = {dimid_obs,dimid_fov,dimid_int_wave};

        // Intensity.
        netcdf_check(nc_l1b,var_radiance_raw = grp_observation.addVar("radiance",ncFloat,dims_spectra_int));
        netcdf_check(nc_l1b,var_radiance_raw.putAtt("long_name","radiance"));
        netcdf_check(nc_l1b,var_radiance_raw.putAtt("units","ph nm-1 s-1 sr-1 m-2"));
        std::vector<std::size_t> chunksize_intensity { 1, ckd->dim_fov, ckd->dim_detector_spec };
        netcdf_check(nc_l1b,var_radiance_raw.setChunking(NcVar::nc_CHUNKED, chunksize_intensity));
        netcdf_check(nc_l1b,var_radiance_raw.setCompression(true, true, 1));

        // Noise on intensity.
        netcdf_check(nc_l1b,var_radiance_raw_noise = grp_observation.addVar("radiance_noise",ncFloat,dims_spectra_int));
        netcdf_check(nc_l1b,var_radiance_raw_noise.putAtt("long_name","noise of I"));
        netcdf_check(nc_l1b,var_radiance_raw_noise.putAtt("units","ph nm-1 s-1 sr-1 m-2"));
        netcdf_check(nc_l1b,var_radiance_raw_noise.setChunking(NcVar::nc_CHUNKED, chunksize_intensity));
        netcdf_check(nc_l1b,var_radiance_raw_noise.setCompression(true, true, 1));

        // Radiance mask
        var_radiance_mask =
          grp_observation.addVar("radiance_mask", ncUbyte, dims_spectra_int);
        var_radiance_mask.putAtt("long_name", "radiance mask");
        constexpr std::array<uint8_t, 2> mask_values { 0, 1 };
        var_radiance_mask.putAtt("flag_values",
                                 netCDF::ncUbyte,
                                 mask_values.size(),
                                 mask_values.data());
        var_radiance_mask.putAtt("flag_meanings", "good bad");
        var_radiance_mask.setChunking(NcVar::nc_CHUNKED, chunksize_intensity);
        var_radiance_mask.setCompression(true, true, 1);
    }

    // NetCDF variables for geolocation.
    if (set->geolocation) {
        dim_utc = 0;

        // Construct the planet.
        earth = make_unique<Planet>(this,set->semi_major_axis,set->semi_minor_axis,set->latitude_tol);
        if (set->demfile.compare("") != 0) {
            handle(earth->include_dem(set->demfile));
        }

        // Create geolocation output.
        vector<NcDim> dims_geo = {dimid_obs,dimid_fov};

        // Latitude.
        netcdf_check(nc_l1b,var_lat = grp_geolocation.addVar("latitude",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_lat.putAtt("standard_name","latitude"));
        netcdf_check(nc_l1b,var_lat.putAtt("long_name","latitude"));
        netcdf_check(nc_l1b,var_lat.putAtt("valid_min",ncFloat,-90.0));
        netcdf_check(nc_l1b,var_lat.putAtt("valid_max",ncFloat,90.0));
        netcdf_check(nc_l1b,var_lat.putAtt("units","degrees_north"));

        // Longitude.
        netcdf_check(nc_l1b,var_lon = grp_geolocation.addVar("longitude",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_lon.putAtt("standard_name","longitude"));
        netcdf_check(nc_l1b,var_lon.putAtt("long_name","longitude"));
        netcdf_check(nc_l1b,var_lon.putAtt("valid_min",ncFloat,-180.0));
        netcdf_check(nc_l1b,var_lon.putAtt("valid_max",ncFloat,180.0));
        netcdf_check(nc_l1b,var_lon.putAtt("units","degrees_east"));

        // Altitude.
        netcdf_check(nc_l1b,var_alt = grp_geolocation.addVar("altitude",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_alt.putAtt("standard_name","altitude"));
        netcdf_check(nc_l1b,var_alt.putAtt("long_name","height above mean sea level"));
        netcdf_check(nc_l1b,var_alt.putAtt("units","m"));
        netcdf_check(nc_l1b,var_alt.putAtt("positive","up"));
        netcdf_check(nc_l1b,var_alt.putAtt("axis","Z"));

        // Viewing zenith angle.
        netcdf_check(nc_l1b,var_vza = grp_geolocation.addVar("sensor_zenith",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_vza.putAtt("long_name","sensor zenith angle"));
        netcdf_check(nc_l1b,var_vza.putAtt("units","degrees"));

        // Viewing azimuth angle.
        netcdf_check(nc_l1b,var_vaa = grp_geolocation.addVar("sensor_azimuth",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_vaa.putAtt("long_name","sensor azimuth angle"));
        netcdf_check(nc_l1b,var_vaa.putAtt("units","degrees"));

        // Solar zenith angle.
        netcdf_check(nc_l1b,var_sza = grp_geolocation.addVar("solar_zenith",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_sza.putAtt("long_name","solar zenith angle"));
        netcdf_check(nc_l1b,var_sza.putAtt("units","degrees"));

        // Solar azimuth angle.
        netcdf_check(nc_l1b,var_saa = grp_geolocation.addVar("solar_azimuth",ncFloat,dims_geo));
        netcdf_check(nc_l1b,var_saa.putAtt("long_name","solar azimuth angle"));
        netcdf_check(nc_l1b,var_saa.putAtt("units","degrees"));

        // Detailed output.
        if (set->l1b_output_level >= 1) {
            netcdf_check(nc_l1b,grp_detailed_output_geolocation = nc_l1b->ncid->addGroup("detailed_output_geolocation"));
            NcDim dimid_vec;
            NcDim dimid_quat;
            netcdf_check(nc_l1b,dimid_vec = grp_detailed_output_geolocation.addDim("vec",DIM_VEC));
            netcdf_check(nc_l1b,dimid_quat = grp_detailed_output_geolocation.addDim("quat",DIM_QUAT));
            vector<NcDim> dims_geo_vec = {dimid_obs,dimid_vec};
            vector<NcDim> dims_geo_quat = {dimid_obs,dimid_quat};
            vector<NcDim> dims_geo_point = {dimid_obs,dimid_fov,dimid_vec};
            netcdf_check(nc_l1b,var_satpos_j2000 = grp_detailed_output_geolocation.addVar("satpos_j2000",ncDouble,dims_geo_vec));
            netcdf_check(nc_l1b,var_rotation_quaternion = grp_detailed_output_geolocation.addVar("rotation_quaternion",ncDouble,dims_geo_quat));
            netcdf_check(nc_l1b,var_t = grp_detailed_output_geolocation.addVar("t",ncDouble,dimid_obs));
            netcdf_check(nc_l1b,var_greenwich_hour_angle = grp_detailed_output_geolocation.addVar("greenwich_hour_angle",ncDouble,dimid_obs));
            netcdf_check(nc_l1b,var_satpos_ecr = grp_detailed_output_geolocation.addVar("satpos_ecr",ncDouble,dims_geo_vec));
            netcdf_check(nc_l1b,var_sunvec_ecr = grp_detailed_output_geolocation.addVar("sunvec_ecr",ncDouble,dims_geo_vec));
            netcdf_check(nc_l1b,var_pointings_ecr = grp_detailed_output_geolocation.addVar("pointings_ecr",ncDouble,dims_geo_point));
            netcdf_check(nc_l1b,var_geoloc_ecr = grp_detailed_output_geolocation.addVar("geoloc_ecr",ncDouble,dims_geo_point));
        }
    }

    // Copy navigation data from the geometry file
    const netCDF::NcFile nc_geom { set->geometry_file, netCDF::NcFile::read };
    NetCDF_object::copyVar(nc_geom.getVar("sza"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });
    NetCDF_object::copyVar(nc_geom.getVar("saa"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });
    NetCDF_object::copyVar(nc_geom.getVar("vza"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });
    NetCDF_object::copyVar(nc_geom.getVar("vaa"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });
    NetCDF_object::copyVar(nc_geom.getVar("lat"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });
    NetCDF_object::copyVar(nc_geom.getVar("lon"),
                           grp_geolocation,
                           { dimid_obs, dimid_fov });

    // Each image is in an individual batch.
    handle(batch_one());

    // These variables are targeted for MPI parallelization. Edit this
    // list if necessary.
    mpi_nc_vars = { { MpiNcVarId::time, var_time } };
    if (ckd->lev > LEVEL_RADCAL) {
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::radiance_raw, var_radiance_raw });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::radiance_raw_noise, var_radiance_raw_noise });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::radiance_mask, var_radiance_mask });
    }
    if (set->geolocation) {
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::lat, var_lat });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::lon, var_lon });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::alt, var_alt });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::vza, var_vza });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::vaa, var_vaa });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::sza, var_sza });
        mpi_nc_vars.emplace_back(std::pair{ MpiNcVarId::saa, var_saa });
    }
    mpi_buffer_starts.resize(mpi_nc_vars.size() + 1, 0);
    mpi_nc_starts.resize(mpi_nc_vars.size());
    mpi_nc_counts.resize(mpi_nc_vars.size());
    mpi_buffers.resize(nl1a_total);
    mpi_buffer_markers.resize(num_procs);
    mpi_requests.resize(nl1a_total, MPI_REQUEST_NULL);
    mpi_buf_multiplicities.resize(nl1a_total, 0);
    mpi_io_timer.reset();
    mpi_final_io_timer.reset();
    // Get the dimensions and size of each netCDF variable in mpi_nc_vars.
    for (int i_var {}; i_var < static_cast<int>(mpi_nc_vars.size()); ++i_var) {
        const int n_dims { mpi_nc_vars.at(i_var).second.getDimCount() };
        mpi_nc_starts.at(i_var).resize(n_dims, 0);
        mpi_nc_counts.at(i_var).resize(n_dims, 1);
        int var_size { 1 };
        for (int i_dim { 1 }; i_dim < n_dims; ++i_dim) {
            mpi_nc_counts.at(i_var).at(i_dim) =
              mpi_nc_vars.at(i_var).second.getDim(i_dim).getSize();
            var_size *= mpi_nc_counts.at(i_var).at(i_dim);
        }
        // Double-precision variables take twice as much space in the buffer.
        var_size *=
          mpi_nc_vars.at(i_var).second.getType().getSize() / sizeof(float);
        mpi_buffer_starts.at(i_var + 1) =
          mpi_buffer_starts.at(i_var) + var_size;
    }
    mpi_buffer_size = mpi_buffer_starts.back();
    mpi_send_size = set->mpi_send_size;
    for (int i_rank { 1 }; i_rank < num_procs; ++i_rank) {
        for (int i_img { mpi_ranges.at(i_rank) };
             i_img < mpi_ranges.at(i_rank + 1);
             i_img += mpi_send_size) {
            // Set every Mth value equal to M where
            // M=mpi_send_size. Be careful at process boundaries (see
            // the comment at mpi_buf_multiplicities).
            mpi_buf_multiplicities.at(i_img) =
              std::min(mpi_send_size, mpi_ranges.at(i_rank + 1) - i_img);
        }
    }
    // The user might have entered a large value for this
    // variable. Reduced it to the true maximum value here (this only
    // matters for output).
    mpi_send_size = std::min(mpi_send_size,
                             *std::max_element(mpi_buf_multiplicities.cbegin(),
                                               mpi_buf_multiplicities.cend()));
    // Initialize mpi_buffer_markers and send out the first
    // asynchronous MPI calls. Most MPI_Irecv calls are initiated at
    // the end of process_batch which is called multiple times.
    if (my_rank == 0) {
        for (int i_rank { 1 }; i_rank < num_procs; ++i_rank) {
            mpi_buffer_markers.at(i_rank) = mpi_ranges.at(i_rank);
        }
        static constexpr int n_buffer_rounds { 2 };
        for (int i {}; i < n_buffer_rounds; ++i) {
            for (int i_rank { 1 }; i_rank < num_procs; ++i_rank) {
                if (mpi_buffer_markers[i_rank] >= mpi_ranges[i_rank + 1]) {
                    break;
                }
                const int first_idx { mpi_buffer_markers.at(i_rank) };
                mpi_buffer_markers.at(i_rank) =
                  first_idx + mpi_buf_multiplicities.at(first_idx);
                mpi_buffers.at(first_idx).resize(
                  mpi_buf_multiplicities.at(first_idx) * mpi_buffer_size);
                MPI_Irecv(mpi_buffers[first_idx].data(),
                          mpi_buf_multiplicities[first_idx] * mpi_buffer_size,
                          MPI_FLOAT,
                          i_rank,
                          first_idx,
                          MPI_COMM_WORLD,
                          &mpi_requests[first_idx]);
            }
        }
    }
    // Buffer size per image in GB
    const double buf_in_GB { static_cast<double>(mpi_buffer_size)
                             * sizeof(float) / 1024 / 1024 / 1024 };
    // Typically, the MPI buffer size on each process is proportional
    // to two times mpi_send_size. Hence the word "Estimated" and the
    // factor of 2.
    writelog(log_info,
             "Estimated MPI buffer size  per process: %.3f GB",
             2 * mpi_send_size * buf_in_GB);
    // Again we say "Estimated" here because the first process will
    // probably hold a slightly higher number of MPI buffers in memory.
    writelog(log_info,
             "Estimated MPI buffer size on process 0: %.3f GB (%.1f%% of the "
             "orbit)",
             (num_procs - 1) * mpi_send_size * buf_in_GB,
             1e2 * (num_procs - 1) * mpi_send_size / nl1a_total);
    if (my_rank == 0) {
        tmp_buf_dbl.resize(mpi_send_size * mpi_buffer_size);
        tmp_buf_float.resize(mpi_send_size * mpi_buffer_size);
    }
    return 0;
} // }}}

// Process one image.
int L1B::process_batch(size_t ibatch, const Calibration_options& opt)
{
    // 1. Interpret index and observation according to L1B batch handling.

    // Each batch is one image.

    // The L1B processor works with one-sized batches, so the observation index is
    // equal to the batch index.
    size_t iobs = ibatch; // More convenient naming, especially if you know SPEXairborne or some biology.
    L1A *l1a = l1a_instances[0]; // Reference to the L1A image.
    // Get the right CKD.
    CKD *ckd = l1a->bin->binned_ckd;

    // All output to be written in the NetCDF.
    vector<double> lat;
    vector<double> lon;
    vector<double> alt;
    vector<double> vza;
    vector<double> vaa;
    vector<double> sza;
    vector<double> saa;
    vector<double> radiance_raw;
    vector<double> radiance_raw_noise;
    vector<bool> radiance_mask;
    vector<double> intensity_radiance;
    vector<double> intensity_radiance_noise;
    vector<double> intens;
    vector<double> intens_noise;

    // 2. Construct the member arrays for output.
    // S-min and S-plus are used inside the L1A module, so are eventuallyh
    // used for intermediate L1X output, so they are allocated when they
    // start making sense (after FOV calibration.
    if (set->geolocation) {
        lat.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        lon.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        alt.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        vza.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        vaa.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        sza.resize(ckd->dim_fov,NC_FILL_DOUBLE);
        saa.resize(ckd->dim_fov,NC_FILL_DOUBLE);
    }
    if (ckd->lev > LEVEL_RADCAL) {
        radiance_raw.resize(ckd->dim_fov_spec_total,NC_FILL_DOUBLE);
        radiance_raw_noise.resize(ckd->dim_fov_spec_total,NC_FILL_DOUBLE);
        radiance_mask.resize(ckd->dim_fov_spec_total, 0);
        intensity_radiance.resize(ckd->dim_fov*dim_int_wave,NC_FILL_DOUBLE);
        intensity_radiance_noise.resize(ckd->dim_fov*dim_int_wave,NC_FILL_DOUBLE);
    }
    // 3. Perform the geo-location.
    if (set->geolocation) {
        const double eR = 6371.0;
        const long double pi = 3.141592653589793;
        const double velocity_lon = l1a->nav_velocity_lon;
        const double velocity_lat = l1a->nav_velocity_lat;
        const double longitude = l1a->nav_longitude;
        const double latitude = l1a->nav_latitude;
        const double h_alt = l1a->nav_altitude;
        const size_t n_act = ckd->dim_fov;
        const double roll_ang = l1a->nav_roll_ang;
        const double h = l1a->nav_altitude;

        std::vector<double> act;
        act = ckd->fov_act_angles;
        for(size_t k = 0; k < n_act; k++){
            lat[k + ibatch * n_act] = latitude + l1a->time * velocity_lat;
            lon[k + ibatch * n_act] = longitude + (360.0/(2 * pi)) * AngToX(roll_ang + act[k], h, eR) / eR;
            vza[k + ibatch * n_act] = (act[k] + roll_ang) +
            (longitude + (360.0/(2 * pi)) * AngToX(roll_ang + act[k], h, eR) / eR); //180-(180-alf-fi)
            sza[k + ibatch * n_act] = 50.0; // The Solar zenith angle, it is considered as 50 at the simple GM model
        }
    }

    // Running pointers. Null-pointers if not used.
    double *radiance_raw_cur = radiance_raw.data();
    //I
    double *intensity_radiance_cur = intensity_radiance.data();
    //I_noise
    double *intensity_radiance_noise_cur = intensity_radiance_noise.data();
    size_t ivp = 0;
    for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
        // This is one spectrum.
        Spectra specs;
        handle(l1a->extract(ifov, opt, specs));
        size_t &nbin = specs.dim;
        // Now, the demodulation starts.
        if (ckd->lev > LEVEL_RADCAL) {
            for (size_t ibin=0 ; ibin<nbin ; ibin++) {
                if (!specs.mask[ibin]) {
                    radiance_raw_cur[ibin] = specs.signal[ibin];
                    intensity_radiance_cur[ibin] = specs.signal[ibin];
                    intensity_radiance_noise_cur[ibin] = specs.noise[ibin];
                } else {
                    radiance_mask[ifov * nbin + ibin] = specs.mask[ibin];
                }
            }
        }
        // Move output pointers.
        if (ckd->lev > LEVEL_RADCAL) {
            // Intensity product and input.
            radiance_raw_cur += ckd->fov_dims_spec[ifov];
            intensity_radiance_cur += dim_int_wave;
            intensity_radiance_noise_cur += dim_int_wave;
        }
    }

    // 5. Write it into the output.
    if (my_rank == 0) {
        vector<size_t> strt_img = {iobs};
        vector<size_t> cnt_img = {1};
        netcdf_check(nc_l1b,var_time.putVar(strt_img,cnt_img,&l1a->time));

        if (ckd->lev > LEVEL_RADCAL) {
            vector<size_t> strt_spec = {iobs,0,0};
            vector<size_t> cnt_spec = {1,ckd->dim_fov,dim_int_wave};
            const auto reverse { [] (const int n_fov,
                                     const int n_wave,
                                     std::vector<double>& data) {
                for (int i {}; i < n_fov; ++i) {
                    for (int j {}; j < n_wave / 2; ++j) {
                        std::swap(data[i * n_wave + j],
                                  data[i * n_wave + n_wave - 1 - j]);
                    }
                }
            } };
            if (ckd->wave_target.front() > ckd->wave_target.back()) {
                reverse(ckd->dim_fov, dim_int_wave, intensity_radiance);
                reverse(ckd->dim_fov, dim_int_wave, intensity_radiance_noise);
            }
            netcdf_check(nc_l1b,var_radiance_raw.putVar(strt_spec,cnt_spec,intensity_radiance.data()));
            netcdf_check(nc_l1b,var_radiance_raw_noise.putVar(strt_spec,cnt_spec,intensity_radiance_noise.data()));
            std::vector<uint8_t> pixel_mask_u8(radiance_mask.size());
            for (size_t i {}; i < pixel_mask_u8.size(); ++i) {
                pixel_mask_u8[i] = static_cast<uint8_t>(radiance_mask[i]);
            }
            netcdf_check(nc_l1b,
                         var_radiance_mask.putVar(strt_spec,
                                                  cnt_spec,
                                                  pixel_mask_u8.data()));
        }
        if (set->geolocation) {
            vector<size_t> cnt_geo = {1,ckd->dim_fov};
            vector<size_t> strt_geo = {iobs,0};
            netcdf_check(nc_l1b,var_lat.putVar(strt_geo,cnt_geo,lat.data()));
            netcdf_check(nc_l1b,var_lon.putVar(strt_geo,cnt_geo,lon.data()));
            netcdf_check(nc_l1b,var_alt.putVar(strt_geo,cnt_geo,alt.data()));
            netcdf_check(nc_l1b,var_vza.putVar(strt_geo,cnt_geo,vza.data()));
            netcdf_check(nc_l1b,var_vaa.putVar(strt_geo,cnt_geo,vaa.data()));
            netcdf_check(nc_l1b,var_sza.putVar(strt_geo,cnt_geo,sza.data()));
            netcdf_check(nc_l1b,var_saa.putVar(strt_geo,cnt_geo,saa.data()));
        }
        netcdf_check(nc_l1b,nc_l1b->ncid->sync());
    }

    if (my_rank > 0) {
        // Start the send process to process 0. Note that
        // mpi_buf_multiplicities is often zero so a buffer is
        // actually only allocated at selected iterations.
        mpi_buffers[iobs].resize(mpi_buf_multiplicities[iobs]
                                 * mpi_buffer_size);
        // The current and maximum allocation size is recorded only
        // for output purposes.
        cur_buf_size += mpi_buffers[iobs].size();
        max_buf_size = std::max(max_buf_size, cur_buf_size);
        // Distance between this iteration and the last iteration when
        // a buffer was allocated. If mpi_buf_multiplicities is zero
        // for this and the previous round but was non-zero two
        // iterations ago then i_shift=2.
        const int i_shift { (static_cast<int>(iobs) - mpi_ranges[my_rank])
                            % mpi_send_size };
        // Index of the last iteration when a buffer was allocated. If
        // a buffer was allocated on this round,
        // i.e. mpi_buf_multiplicities is non-zero for this round then
        // i_active_idx=iobs.
        const int i_active_idx { static_cast<int>(iobs) - i_shift };
        // Where to start writing targeted variables (DoLP, AoLP, ...)
        // in the MPI buffer. There might already be data from
        // previous iterations so need to shift the starting index
        // accordingly.
        const int buf_start_idx { i_shift * mpi_buffer_size };
        auto buf_start { mpi_buffers[i_active_idx].begin() + buf_start_idx };
        // Main loop where we copy each netCDF variable into the MPI
        // buffer one by one.
        for (int i_var {}; i_var < static_cast<int>(mpi_nc_vars.size());
             ++i_var) {
            const auto copy { [this, &buf_start, i_var](const auto& source) {
                std::transform(source.cbegin(),
                               source.cend(),
                               buf_start + mpi_buffer_starts[i_var],
                               [](const double x) {
                                   return static_cast<float>(x);
                               });
            } };
            switch (mpi_nc_vars[i_var].first) {
            case MpiNcVarId::time:
                // We are copying a double precision variable into a
                // single precision buffer.
                std::memcpy(&mpi_buffers[i_active_idx][
                              buf_start_idx + mpi_buffer_starts[i_var]],
                            &l1a->time,
                            sizeof(double));
                break;
            case MpiNcVarId::radiance_raw:
                copy(intensity_radiance);
                break;
            case MpiNcVarId::radiance_raw_noise:
                copy(intensity_radiance_noise);
                break;
            case MpiNcVarId::radiance_mask:
                copy(radiance_mask);
                break;
            case MpiNcVarId::intens:
                copy(intens);
                break;
            case MpiNcVarId::intens_noise:
                copy(intens_noise);
                break;
            case MpiNcVarId::lat:
                copy(lat);
                break;
            case MpiNcVarId::lon:
                copy(lon);
                break;
            case MpiNcVarId::alt:
                copy(alt);
                break;
            case MpiNcVarId::vza:
                copy(vza);
                break;
            case MpiNcVarId::vaa:
                copy(vaa);
                break;
            case MpiNcVarId::sza:
                copy(sza);
                break;
            case MpiNcVarId::saa:
                copy(saa);
                break;
            }
        }
        // If this is the final iteration for the given MPI process or
        // for the next iteration mpi_buf_multiplicities is non-zero
        // again, it's time to send the buffer contents to process 0.
        if (iobs == mpi_ranges[my_rank + 1] - 1
            || mpi_buf_multiplicities[iobs + 1] > 0) {
            MPI_Isend(mpi_buffers[i_active_idx].data(),
                      mpi_buf_multiplicities[i_active_idx] * mpi_buffer_size,
                      MPI_FLOAT,
                      0,
                      i_active_idx,
                      MPI_COMM_WORLD,
                      &mpi_requests[i_active_idx]);
        }
        // Check if process 0 has recieved any of the previous
        // requests. If so, deallocate those MPI buffers.
        for (int i_img { mpi_ranges[my_rank] };
             i_img < mpi_ranges[my_rank + 1];
             i_img += mpi_send_size) {
            if (mpi_requests[i_img] != MPI_REQUEST_NULL) {
                MPI_Status status {};
                int flag {};
                MPI_Test(&mpi_requests[i_img], &flag, &status);
                if (flag) {
                    cur_buf_size -= mpi_buffers[i_img].size();
                    mpi_buffers[i_img] = std::vector<float>{};
                }
            }
        }
    } else {
        mpi_io_timer.start();
        const bool last_round { iobs == mpi_ranges[1] - 1 };
        if (last_round) {
            mpi_final_io_timer.start();
        }
        // Tests if any data was received in the main double loop.
        bool active_run {};
        do {
            active_run = false;
            for (int i_rank { 1 }; i_rank < num_procs; ++i_rank) {
                for (int i_img { mpi_ranges[i_rank] };
                     i_img < mpi_ranges[i_rank + 1];
                     i_img += mpi_send_size) {
                    // MPI_REQUEST_NULL signals that data has already
                    // been received from that buffer. Finish this
                    // iteration.
                    if (mpi_requests[i_img] == MPI_REQUEST_NULL) {
                        continue;
                    }
                    active_run = true;
                    MPI_Status status;
                    int flag {};
                    // Test if the MPI_Isend operation for image i_img
                    // has been completed.
                    MPI_Test(&mpi_requests[i_img], &flag, &status);
                    // If not, try again later but finish this
                    // iteration for now.
                    if (!flag) {
                        continue;
                    }
                    // MPI buffer i_img contains data. Write that data
                    // to output file. First treat image_time which is
                    // the only double-precision variable.
                    for (int i_count {};
                         i_count < mpi_buf_multiplicities[i_img];
                         ++i_count) {
                        // i_count iterates over the number of images
                        // received in this batch.
                        const int start_buf { i_count * mpi_buffer_size
                                              + mpi_buffer_starts.front() };
                        std::memcpy(&tmp_buf_dbl[i_count],
                                    &mpi_buffers[i_img][start_buf],
                                    sizeof(double));
                    }
                    mpi_nc_starts.front().front() = i_img;
                    mpi_nc_counts.front().front() =
                      mpi_buf_multiplicities[i_img];
                    mpi_nc_vars.front().second.putVar(
                      mpi_nc_starts.front(),
                      mpi_nc_counts.front(),
                      tmp_buf_dbl.data());
                    // Next treat all the other variables, which are
                    // single-precision.
                    for (int i_var { 1 };
                         i_var < static_cast<int>(mpi_nc_vars.size());
                         ++i_var) {
                        for (int i_count {};
                             i_count < mpi_buf_multiplicities[i_img];
                             ++i_count) {
                            const int var_size { mpi_buffer_starts[i_var + 1]
                                                 - mpi_buffer_starts[i_var] };
                            const int start_tmp { i_count * var_size };
                            const int start_buf { i_count * mpi_buffer_size
                                                  + mpi_buffer_starts[i_var] };
                            for (int i {}; i < var_size; ++i) {
                                tmp_buf_float[i + start_tmp] =
                                  mpi_buffers[i_img][i + start_buf];
                            }
                        }
                        mpi_nc_starts[i_var].front() =
                          i_img;
                        mpi_nc_counts[i_var].front() =
                          mpi_buf_multiplicities[i_img];
                        mpi_nc_vars[i_var].second.putVar(
                          mpi_nc_starts[i_var],
                          mpi_nc_counts[i_var],
                          tmp_buf_float.data());
                    }
                    // Mark operation as complete (see the conditional
                    // at the beginning of this loop).
                    mpi_requests[i_img] = MPI_REQUEST_NULL;
                    // Free the used up MPI buffer and continue to
                    // allocate a new buffer and send out a
                    // corresponding asynchronous receive request.
                    mpi_buffers[i_img] = std::vector<float>{};
                    // If the buffer marker has passed the image range
                    // assigned to a given process then all buffers
                    // have already been allocated.
                    if (mpi_buffer_markers[i_rank] >= mpi_ranges[i_rank + 1]) {
                        continue;
                    }
                    const int next_img { mpi_buffer_markers[i_rank] };
                    mpi_buffers[next_img].resize(mpi_buf_multiplicities[i_img]
                                                 * mpi_buffer_size);
                    MPI_Irecv(mpi_buffers[next_img].data(),
                              mpi_buf_multiplicities[i_img] * mpi_buffer_size,
                              MPI_FLOAT,
                              i_rank,
                              next_img,
                              MPI_COMM_WORLD,
                              &mpi_requests[next_img]);
                    mpi_buffer_markers[i_rank] += mpi_buf_multiplicities[i_img];
                }
            }
            // The double loop runs as many times as necessary until
            // there are no more active MPI requests and all data has
            // been received by process 0.
        } while (last_round && active_run);
        if (last_round) {
            mpi_final_io_timer.stop();
        }
        mpi_io_timer.stop();
    }
    return 0;
} // }}}

// Finalize (write the creation time).
int L1B::process_finalize( // {{{
)
{
    if (num_procs > 1) {
        if (my_rank > 0) {
            const int n_peak_images { static_cast<int>(max_buf_size
                                                       / mpi_buffer_size) };
            MPI_Send(&n_peak_images, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
            // Delete the dummy output files
            nc_l1b->close();
            std::remove(nc_l1b->filename.c_str());
        } else {
            std::vector<int> n_peak_images(num_procs);
            for (int i_rank { 1 }; i_rank < num_procs; ++i_rank) {
                MPI_Recv(&n_peak_images.at(i_rank),
                         1,
                         MPI_INT,
                         i_rank,
                         i_rank,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                writelog(log_info,
                         "Maximum number of images held in the MPI buffer on "
                         "process %3i: %i",
                         i_rank,
                         n_peak_images.at(i_rank));
            }
        }
        writelog(log_info,
                 "MPI communication time on process 0: %.3f s (%i calls total)",
                 mpi_io_timer.totalWallTime(),
                 mpi_io_timer.totalNumberOfCalls());
        writelog(log_info,
                 "Final communication step on process 0: %.3f s",
                 mpi_final_io_timer.totalWallTime());
    }
    return 0;
} // }}}
