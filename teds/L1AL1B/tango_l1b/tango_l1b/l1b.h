// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

namespace tango {

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
}; // }}}

} // namespace tango
