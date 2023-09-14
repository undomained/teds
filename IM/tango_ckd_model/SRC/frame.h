#ifndef FRAME_H
#define FRAME_H

#include "header.h"
#include "logger.h"

class Settings_l1b;
class Settings_isrf;
class Settings_geo;
class UTC;
class CKD;
class L1X_inputfile;
class Planet;

// Class for one frame and everything that is being done with it.
class Frame : public Logger { // {{{

    public:
    Frame(
        Logger *creator
    );
    ~Frame();

    int resample(
        CKD *ckd, // Calibration key data.
        int fct // Enahncement factor.
    );

    int convert_units(
        CKD *ckd
    );

    int read_l1x(
        size_t iframe,
        CKD *ckd,
        L1X_inputfile *l1x_inputfile
    );

    int write_truth_l1b(
        CKD *ckd,
        Settings_l1b *set
    );

    int isrf_interpolate(CKD *ckd, Settings_isrf *set);

    int interpolate_truth(
        CKD *ckd
    );

    int modulate(
        CKD *ckd
    );

    int uncalibrate_spectra(
        CKD *ckd
    );

    int draw_on_detector(
        CKD *ckd
    );

    int apply_straylight(const Settings_main& set, CKD *ckd);

    int apply_prnu(
        CKD *ckd
    );

    int apply_nonlinearity(
        CKD *ckd
    );

    int apply_dark_current(
        CKD *ckd
    );

    int apply_dark_offset(
        CKD *ckd
    );

    // Contents.
    size_t l1x_iframe;

    // Image attributes.
    double image_time; // Interpolation time axis.
    uint32_t seconds; // Absolute time in seconds from 1970.
    int microseconds; // Microseconds (0-999999).
    // No binning table, because we do not support binning tables except at the very end.
    uint16_t nr_coadditions; // Number of co-additions in final L1A file.
    double exposure_time; // Detector exposure time.

    // Engineering data.
    double detector_temperature; // Detector temperature (for the dark).

    // Navigation data, still in array format. Vectors and Quaternions will
    // be created when using them in calculations. Here, they stay in the
    // format so that they can directly be written ito the output NetCDF.
    bool navi_exists;
    vector<double> orb_pos; // Satellite position.
    vector<double> orb_vel; // Satellite velocity.
    vector<double> att_quat; // Attitude quaternion.

    // True image.
    size_t dim_spec_truth;
    size_t dim_spat_truth;
    vector<double> wavelength;
    vector<double> intens;

    // SPEX spectra (calibrated and uncalibrated).
    vector<double> spex_spectra; // S- and S+.

    // Detector image.
    vector<double> image;
    vector<double> image_with_current;

    vector<int> image_ints;

}; // }}}

#endif
