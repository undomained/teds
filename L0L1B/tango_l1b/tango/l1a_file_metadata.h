// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1A_FILE_META_H
#define L1A_FILE_META_H

#include "header.h"
#include "logger.h"

// Forward declaration.
class Vector;
class Quaternion;
class NetCDF_object;

// L1A file information.
// This is information that is specific for one L1A file, which need not be exactly
// coupled to one image. For instance, engineering data has to be interpolated for each image.
// If multiple images exist in one file (such as in flight L1A), we would not like to repeat
// reading the engineering data for each image.
class L1A_file_metadata : public Logger { // {{{

    public:
    // Constructor.
    L1A_file_metadata(
        Logger *creator
    );
    shared_ptr<NetCDF_object> nc; // NetCDF object of the L1A file.
    // After processing the metadata, and instances of L1A_file_metadata no longer exist,
    // the ownership of the NetCDF object will be transferred to the L1A structure.
    string filename; // Why not? Of course this is set->l1a_files[ifile], but now, the
    // L1A_file_metadata object can work on its own.
    size_t file_id; // For possible recognition, relevant for the frames L1A manager if the processor wants to know what belongs to what.
    size_t nframe; // Number of frames in the file, so important for organizational issues, that we save this number as soon as possible.
    size_t npix; // Number of pixels for image (maximum).

    // Image time attributes (to be copied to L1B).
    // This is exactly the same as the orbit metadata structure.
    // Only here, it is just added to the members of L1A file metadata.
    // There is less notion of rubbish here, because this class is a
    // garbage dump anyway, and much more hidden than the processor or
    // the L1A manager.
    string startDirection; // Ascending or descending.
    string endDirection; // Ascending or descending.
    long long orbit_number;
    string time_coverage_start;
    string time_coverage_stop;
    string time_reference;

    // Metadata with time axis different from the images.

    // Engineering data.
    size_t eng_ntime;
    vector<double> eng_time;
    vector<double> eng_temp_detector;
    // We do not (yet) use the temp_optics.

    // Navigation data (if requested).
    // For some reason, the navigation data has two time axes, altough it is hardcoded that they
    // have the same number of elements.
    bool l1a_navigation; // Request flag for navigation data.
    size_t nav_ntime;
    vector<double> nav_att_time;
    vector<Quaternion> nav_att_quat;
    vector<double> nav_orb_time;
    vector<Vector> nav_orb_pos;
    vector<Vector> nav_orb_vel;
    
    // An added challenge that these do not interpolate straight forward. The position and
    // velocity use an interpolation where the velocity is the derivative of the position, with
    // Earth rotation correction. And Quaternions are interpolated symmetrical and naturally
    // normalized (see eventual issues in the appendix in the ATBD.)

    // GSE-data (optional, and totally absent for flight L1A).
    // This group exists per file.
    size_t dim_refspec; // Number of wavelength of a reference spectrum (only for optics measurements).
    double illumination_level; // Scalar illumination level for detector measurements.
    uint8_t viewport; // Illuminated viewports.
    vector<double> refspec_wavelength; // Dimension refspec.
    vector<double> refspec_radiance; // Dimension refspec.
    double refspec_dolp; // Constant over the entire spectrum.
    double refspec_aolp; // Constant over the entire spectrum.
    double act_angle; // Across-track rotation angle of rotation stage.
    double alt_angle; // Along-track rotation angle of rotation stage.

    // Global attribute only for L1X input (developed L1A).
    int il1x_start; // L1X maturity level. Fillvalue for normal L1A.

    // L1X dimension check. There is no CKD when you read the input.
    // Desired dimensions of unbinned images or spectra should be
    // communicated to the L1A and checked with the CKD. After that, the
    // dimensions are no longer read during the read_image routine.
    size_t l1x_dim_detector_spec = NC_FILL_UINT64; // Ensures that any mistake will result in a segmentation fault.
    size_t l1x_dim_detector_spat = NC_FILL_UINT64; // Same idea.
    size_t l1x_dim_fov = NC_FILL_UINT64; // Same idea.

    // Read routine.
    int read_metadata();

}; // }}}

#endif
