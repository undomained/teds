// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1A_H
#define L1A_H

#include "header.h"
#include "logger.h"
#include "vector.h"

// Forward declaration.
class Netcdf_object;
class Calibration_options;
class CKD;
class L1A_file_metadata;
class L1X;
class Binningtable;

// Help struct: Extracted spectra.
struct Spectra {
    size_t dim;
    vector<double> signal;
    vector<double> noise;
    vector<bool> mask;
};

// Class for one image and everything that is being done with it.
class L1A : public Logger { // {{{

    public:
    L1A(
        Logger *creator
    );
    ~L1A();

    L1A_file_metadata *file_meta; // Pointer to file-based metadata. Pointer will expire after initialization.
    shared_ptr<NetCDF_object> nc; // Pointer to NetCDF object that belongs to this L1A image.
    // If each L1A-image has its own file, this pointer is unique. If one file is used for all
    // images, this pointer is shared among all L1A instances and the corresponding L1A manager.
    string filename; // Why not? Of course this is set->l1a_files[ifile], but now, the
    size_t file_id; // For possible recognition, relevant for the frames L1A manager if the processor wants to know what belongs to what.
    size_t frame_offset; // Frame number where the current observation starts (always 0 for calibration L1A).
    size_t nframe_inp; // Number of frames in NetCDF to be coadded.

    // Chunk that is actually read.
    size_t ipix_start = 0;
    size_t ipix_end = 0;

    // Binning table.
    // As metadata, all binning table identifiers (a number) are read. Then, the L1A manager
    // assesses which unique binning tables there are, creates them as local vector and then
    // distributes them. Therefore, the object itself is a shared pointer. The binning table can
    // be large with interpolation instructions.
    uint8_t binning_table_id; // Number that represents the binning table.
    shared_ptr<Binningtable> bin; // Pointer to binning table.

    // From the settings.
    bool l1a_navigation; // Request flag for navigation data.

    //The follows has been added for the simplified GM
    double nav_latitude;
    double nav_longitude;
    double nav_velocity_lat;
    double nav_velocity_lon;
    double nav_altitude;
    double nav_roll_ang;
    //The end of the added for the simple GM

    // Main data that is read.
    vector<double> image_read; // Fraction of the image that is read. Usually the full image, but not always.
    vector<double> noise_read; // Same for the noise, though it is not read but calculated.

    // Pointers such that original pixel index can be used.
    double *image; // Pointer to the image such that the normal pixel index can be used, also for fragmented images.
    double *noise; // It works the same for noise.

    // Temporarily needs two image during noise evaluation.
    vector<double> image_with_dark_current_read; // This will be subtracted with a delay. Keep it zero if dark is skipped.
    double *image_with_dark_current; // Also for noise calibration application during calibrate detector.

    // Pixel mask, always fully defined.
    vector<bool> pixelmask; // Dimensions (dim_detector_spat,dim_detector_spec). Acquired from mask file, flag for false=living, true=dead by L1A mask.

    // Image attributes.
    int index; // Image identifier.
    double time; // Time of the measurement used for interpolation.
    // Absolute time for solar model.
    int seconds; // Number of seconds since Greenwich new year 1970.
    double secondfraction; // Additional fraction of a second to be added (just to add precision to the absolute time.

    // Meta-data, eventuallhy averaged or added when co-adding.
    double exposure_time; // Exposure time (eventually averaged).
    int nr_coadditions; // Total number of co-additions.

    // GSE-data (optional, and totally absent for flight L1A).
    size_t dim_refspec; // Number of wavelength of a reference spectrum (only for optics measurements).
    double illumination_level; // Scalar illumination level for detector measurements.
    uint8_t viewport; // Illuminated viewports.
    vector<double> refspec_wavelength; // Dimension refspec.
    vector<double> refspec_radiance; // Dimension refspec.
    double refspec_dolp; // Constant over the entire spectrum.
    double refspec_aolp; // Constant over the entire spectrum.
    double act_angle; // Across-track rotation angle of rotation stage.
    double alt_angle; // Along-track rotation angle of rotation stage.

    // Naviagion data.
    // The attitude quaternion and satellite position can be interpolated directly.
    Quaternion attitude_quaternion; // Rotation quaternion for attitude.
    Vector satellite_position; // Interpolated satellite position.
    // Both are in J2000 system.
    Vector satellite_velocity; // Only for L1X.

    // L1X output.
    vector<shared_ptr<L1X>> l1x; // L1X instances for each step. Null-pointers for non-written steps.
    size_t l1x_iframe; // Which frame in the L1X is owned by this L1A instance. It is always only one frame because co-adding is the first thing that happens.

    Spectra *l1x_specs;
    Spectra *l1x_specs_target;

    // L1X input.
    int il1x_start;
    // When spectra are read, they are stashed into this array.
    // This is needed because the L1A usually has no spectra as members.
    // The idea is to read the L1X spectra and stash them. When a spectrum
    // is requested, it is retrieved from this stash rather than acquired
    // from the image, which should not exist if you have L1X of FOV
    // level or higher.
    vector<Spectra> l1x_stashed_spectra;
    // L1X dimension check. There is no CKD when you read the input.
    // Desired dimensions of unbinned images or spectra should be
    // communicated to the L1A and checked with the CKD. After that, the
    // dimensions are no longer read during the read_image routine.
    size_t l1x_dim_detector_spec;
    size_t l1x_dim_detector_spat;
    size_t l1x_dim_fov;
    // TODO: Dimensions change when binning procedure has been updated.

    // Read routines.
    int read_metadata();
    int read_image(
        size_t a_ipix_start,
        size_t a_ipix_end
    );
    int remove_image();

    // Calibrated an L1A image for detector issues.
    int calibrate_detector(
        CKD *ckd, // Calibration key data.
        Calibration_options &opt // Detector calibration options.
    );

    // Extracts a spectrum from a detector image.
    // This takes responsibility of the FOV-calibration, but only extracts
    // one spectrum at a time. If radiometric calibration key data is present,
    // it will be radiometrically calibrated.
    int extract(
        size_t ifov, // Field-of-view index.
        const Calibration_options& opt,
        Spectra &specs // Output spectra.
    );

}; // }}}

#endif
