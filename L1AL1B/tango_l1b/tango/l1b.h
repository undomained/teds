// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1B_H
#define L1B_H

#include "header.h"
#include "settings_proc.h"
#include "parallel.h"
#include "processor.h"

// Forward declaration.
class NetCDF_object;
class Planet;
class CKD;

class Settings_l1b : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_l1b(
        Logger *creator
    );
    // Specific settings.
    string outputfile = ""; // Output L1B file.
    size_t order = 0; // B-spline order of wavelength dependence of q and u.
    vector<double> knots; // Wavelengths of knot positions in B-spline.
    double wave_min = NC_FILL_DOUBLE; // Wavelength to start the demodulation window.
    double wave_max = NC_FILL_DOUBLE; // Wavelength to end the demodulation window.
    vector<double> l1b_wavelength; // Central wavelengths of polarimetric product.
    vector<double> intensity_wavelength; // Wavelengths for intensity-only grid.
    double resolving_power = NC_FILL_DOUBLE; // Full-width half-maximum of Gauss convolution of polarimetric L1B product.
    double gauss_range = NC_FILL_DOUBLE; // Range of Gauss in FWHMs on both sides.
    bool demod_noise = true; // Flag to also calculate noise propagation on demodulation products.
    bool &geolocation = l1a_navigation; // Request for geo-location, auto-linked with request for L1A navigation.
    string utcfile = ""; // UTC Time-difference file.
    string demfile = ""; // Detailed elevation map. Leave empty for placeholder DEM function.
    double semi_major_axis = NC_FILL_DOUBLE; // Semi-major axis of the earth.
    double semi_minor_axis = NC_FILL_DOUBLE; // Semi-minor axis of the earth.
    double latitude_tol = NC_FILL_DOUBLE; // Tolerance for iterative calculation of latitude during geolocation (radians).
    double mountain = NC_FILL_DOUBLE; // Highest mountain to expect.
    double movedistance = NC_FILL_DOUBLE; // Distance to move in one iteration avoiding to skip mountains.
    double extremeweight = NC_FILL_DOUBLE; // Clip on weight factor for guess during geolocation convergence.
    double geolocation_tol = NC_FILL_DOUBLE; // Tolerance for convergence along line of sight (length units).
    int l1b_output_level = 0; // Level of detailed output. Own implementation, because output file is not the CKD file.
    // Relative workload on the first MPI process compared to the
    // other processes. Example: if this value is 0.9 (i.e. 90%) then
    // with 10 MPI processes and 1000 images process 0 is assigned the
    // first 90 images. The reason process 0 should be dealt less
    // images is that it also needs to perform I/O operations.
    double first_proc_rel_workload { 0.9 };
    // Number of images each MPI task processes before sending them to
    // process 0. Default is 10, i.e. an asynchronous MPI send is
    // initiated after every 10th image. Increasing this will reduce
    // the amount of communication but increase memory usage on all
    // MPI processes.
    int mpi_send_size { 10 };
    // NetCDF file containing geolocation data (lat, lon, sza, saa, vza, vaa)
    std::string geometry_file {};

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

struct Bin_initialization { // {{{
    uint8_t binning_table_id; // Identifier of the binning table that corresponds to this initialization structure.
    vector<size_t> ibin_start; // Binned spectral index where the retrieval window starts.
    vector<size_t> imeas_start_fov; // Starting index per FOV on the measurement (retrieval window) / FOV folded dimension.
    vector<size_t> nmeas_fov; // Number of measurements in the retrieval window per FOV.
    vector<bool> pol_mask; // Mask arising from holes in CKD (relevant for raw polcal calculations).
    vector<double> jacobians; // Inversion jacobian matrix for each viewport and FOV.
    vector<double> gausses; // All gausses times delta wavelength for convolution.
    vector<size_t> gauss_imeas_start; // Starting indices of the Gauss domain.
    vector<size_t> gauss_imeas_end; // Ending indices of the Gauss domain.
    vector<double> bmat; // B-spline matrices for each viewport and FOV.
    vector<size_t> spline_imeas_start; // Starting indices of the splines.
    vector<size_t> spline_imeas_end; // Ending indices of the splines.
}; // }}}

class L1B : public Processor { // {{{

    // A public constructor.
    public:
    L1B(
        Logger *creator,
        CKD *ckd_arg
    );
    ~L1B();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_l1b> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(size_t ibatch, const Calibration_options& opt) override;
    int process_finalize() override;

    private:

    // An extra dimension.
    size_t dim_pol_wave; // Number of wavelenghts in (polarization) L1B output.
    size_t dim_int_wave; // Number of wavelengths in intensity-only L1B output.

    // Output.
    unique_ptr<NetCDF_object> nc_l1b;

    NcVar var_time;
    // Raw Radiance.
    NcVar var_radiance_raw;
    NcVar var_radiance_raw_noise;
    NcVar var_radiance_mask;
    // Demodulation.
    NcVar var_intens;
    NcVar var_small_q;
    NcVar var_small_u;
    NcVar var_dolp;
    NcVar var_aolp;
    NcVar var_intens_noise;
    NcVar var_small_q_noise;
    NcVar var_small_u_noise;
    NcVar var_dolp_noise;
    NcVar var_aolp_noise;
    // Geolocation.
    NcVar var_lat;
    NcVar var_lon;
    NcVar var_alt;
    // Geometry.
    NcVar var_vza;
    NcVar var_vaa;
    // Solar geometry.
    NcVar var_sza;
    NcVar var_saa;
    // Detailed output.
    NcGroup grp_detailed_output_geolocation;
    NcVar var_satpos_j2000;
    NcVar var_rotation_quaternion;
    NcVar var_t;
    NcVar var_greenwich_hour_angle;
    NcVar var_satpos_ecr;
    NcVar var_sunvec_ecr;
    NcVar var_pointings_ecr;
    NcVar var_geoloc_ecr;

    size_t ntable;
    vector<Bin_initialization> bin_inits; // Binning-table dependent initializations.

    size_t nspline; // Number of spline functions in basis set.
    size_t nstate; // Number of state parameters of each of the inversion.

    // UTC file contents (for geolocation).
    size_t dim_utc;
    vector<int> mjd_utc;
    vector<double> tdiff_utc;

    // The planet.
    unique_ptr<Planet> earth; // Just a silly name of a variable from type Planet.

    // We use these enums to connect netCDF variables with the
    // corresponding data, e.g. var_dolp and dolp.
    enum class MpiNcVarId
    {
        time,
        radiance_raw,
        radiance_raw_noise,
        radiance_mask,
        intens,
        intens_noise,
        small_q,
        small_q_noise,
        small_u,
        small_u_noise,
        dolp,
        dolp_noise,
        aolp,
        aolp_noise,
        lat,
        lon,
        alt,
        vza,
        vaa,
        sza,
        saa,
    };
    // These netCDF variables are targeted for MPI parallelization.
    std::vector<std::pair<MpiNcVarId, netCDF::NcVar>> mpi_nc_vars {};
    // Each MPI buffer in mpi_buffers corresponds to one image,
    // although not all nl1a_total images are allocated at the same
    // time. Each MPI buffer contains the data defined by mpi_nc_vars,
    // i.e. variables such as DoLP, AoLP, etc. are serialized into a
    // single buffer per iteration (image).
    std::vector<std::vector<float>> mpi_buffers {};
    // The starting index of each quantity (DoLP, AoLP, ...) in an MPI buffer
    std::vector<int> mpi_buffer_starts {};
    // When writing data to file, defines the starting indices of each
    // netCDF variable (the startp argument of netCDF::NcVar::putVar)
    std::vector<std::vector<size_t>> mpi_nc_starts {};
    // When writing data to file, defines the number of indices along
    // each dimension for the netCDF variables (the countp argument of
    // netCDF::NcVar::putVar)
    std::vector<std::vector<size_t>> mpi_nc_counts {};
    // Size of an MPI buffer per image. Note that each member of
    // mpi_buffer is resized to an integer multiple (greater than 1)
    // of this variable because it's more efficient to send multiple
    // images in a batch, i.e. not initiate an asynchronous MPI
    // command after each iteration.
    int mpi_buffer_size {};
    // For process 0 (my_rank=0), keeps track of images for which the
    // data has already been received. The total number of those
    // images is nl1a_total minus those assigned to process 0. For all
    // other processes keeps track of MPI_Isend operations that have
    // successfully completed and thus signals that the corresponding
    // memory may be released.
    std::vector<MPI_Request> mpi_requests {};
    // Keeps track of the highest image index for which a buffer has
    // been allocated for the given MPI process. For instance, with
    // 9000 images and 20 MPI processes, let's say process 15 has been
    // assigned images 6750...7200 but only the MPI buffers for images
    // 6750...7014 have been allocated because it's too expensive to
    // keep the whole orbit in memory. Then
    // mpi_buffer_markers[14]=7014. This value is increased as the
    // calculation progresses until all buffers have been allocated
    // and all data transmitted to process 0.
    std::vector<int> mpi_buffer_markers {};
    // Same meaning as Settings_l1b::mpi_send_size
    int mpi_send_size { 1 };
    // Determines the size of the MPI buffer allocated on a given
    // iteration inside the main loop over images. Its values are
    // typically either mpi_send_size or 0. In the following example
    // the columns are R - MPI rank, N - image index, M -
    // mpi_buf_multiplicities[N].
    //
    // R  N  M
    // ...
    // 1 34  0
    // 2 35  3 # Allocate MPI_buffer with size 3 * mpi_buffer_size
    // 2 36  0 # Do not allocate new buffer; write data into previously
    // 2 37  0 #    allocated buffer
    // 2 38  2 # Allocate MPI_buffer with size 2 * mpi_buffer_size and not
    // 2 39  0 #    2 * mpi_buffer_size because there aren't enough images left.
    // 3 40  3
    // ...
    //
    // For image 35 an MPI buffer is allocated with size 3 *
    // mpi_buffer_size (meaning mpi_send_size=3). In subsequent
    // iterations a new buffer is not allocated and instead data is
    // written into the previously allocated buffer. At the end of
    // iteration 37 data is sent to process 0 and soon after the
    // buffer is released. At iteration 38 a new buffer is allocated
    // but with the size 2 * mpi_buffer_size because there are only
    // two images left (image #40 is assigned to the next MPI
    // process).
    std::vector<int> mpi_buf_multiplicities {};
    // Keeps track of the time spent by process 0 collecting data from
    // other processes.
    gadfit::Timer mpi_io_timer {};
    // Keeps track of the time spent by process 0 in the final step
    // collecting data from other processes.
    gadfit::Timer mpi_final_io_timer {};
    // Current number of MPI buffers allocated on this MPI process
    size_t cur_buf_size {};
    // Maximum number of MPI buffers allocated on this MPI process
    size_t max_buf_size {};
    // In mpi_ranges, data corresponding to different images is stored
    // sequentially. This means that all data corresponding to a
    // single variable such as dolp is stored non-contiguously. We use
    // tmp_buf_dbl/float for rearranging the data layout before
    // writing it to file.
    std::vector<double> tmp_buf_dbl {};
    std::vector<float> tmp_buf_float {};
}; // }}}

#endif
