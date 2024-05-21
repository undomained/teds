#ifndef SETTINGS_L1B_H
#define SETTINGS_L1B_H

#include "header.h"
#include "settings.h"

class NetCDF_object;
class CKD;

// L1B settings.
class Settings_l1b : public Settings {

    public:

    // Define settings here.
    bool execute = false; // Flag for writing truth L1B file.
    string outputfile; // Output L1B file.
    vector<double> l1b_wavelength; // Central wavelengths of polarimetric product.
    vector<double> intensity_wavelength; // Wavelengths for intensity-only grid.
    double resolving_power; // Full-width half-maximum of Gauss convolution of polarimetric L1B product.
    double gauss_range; // Range of Gauss in FWHMs on both sides.

    // Constructor.
    Settings_l1b(
        Logger *creator
    );
    ~Settings_l1b(); // Destructor.

    // L1B NetCDF stuff.
    unique_ptr<NetCDF_object> nc_l1b;
    NcVar var_image_time;
    NcVar var_radiance_raw;
    NcVar var_intens;
    NcVar var_small_q;
    NcVar var_small_u;
    NcVar var_dolp;
    NcVar var_aolp;

    size_t dim_int_wave; // Number of intensity wavelengths.
    size_t dim_pol_wave; // Number of polarization wavelengths.

    // Overwritten common settings.
    protected:
    int init_common(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) override;

};

#endif
