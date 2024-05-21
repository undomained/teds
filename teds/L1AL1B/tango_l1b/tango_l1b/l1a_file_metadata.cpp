// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "logger.h"
#include "vector.h"
#include "netcdf_object.h"
#include "l1a_file_metadata.h"

#include <algorithm>

namespace tango {

L1A_file_metadata::L1A_file_metadata( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}

// Reads metadata of one file.
// This routine does not know what image is interesting, but reads important file-specific
// information, such as the number of frames, as well as the metadata with different time axes.
int L1A_file_metadata::read_metadata( // {{{
)
{
    // Make a shared pointer to the NetCDF object, will be transferred to the L1A images after
    // L1A manager initialization.
    nc = make_shared<NetCDF_object>(this);
    // Open the file.
    handle(nc->open(filename,NcFile::read));

    { // Search global attributes for orbit data.
        multimap<string,NcGroupAtt> mp;
        netcdf_check(nc,mp = nc->ncid->getAtts());
        multimap<string,NcGroupAtt>::iterator search;
        search = mp.find("startDirection");
        if (search == mp.end()) startDirection = "Unknown";
        else {
            netcdf_check(nc,search->second.getValues(startDirection));
        }
        search = mp.find("endDirection");
        if (search == mp.end()) endDirection = "Unknown";
        else {
            netcdf_check(nc,search->second.getValues(endDirection));
        }
        search = mp.find("orbit_number");
        if (search == mp.end()) orbit_number = NC_FILL_INT64;
        else {
            netcdf_check(nc,search->second.getValues(&orbit_number));
        }
        search = mp.find("time_coverage_start");
        if (search == mp.end()) time_coverage_start = "Unknown";
        else {
            netcdf_check(nc,search->second.getValues(time_coverage_start));
        }
        search = mp.find("time_coverage_stop");
        if (search == mp.end()) time_coverage_stop = "Unknown";
        else {
            netcdf_check(nc,search->second.getValues(time_coverage_stop));
        }
        // L1X maturity level. Does not exist for L1A, but does for
        // L1X input files.
        search = mp.find("l1x_maturity_level");
        if (search == mp.end()) il1x_start = L1X_L1A;
        else {
            netcdf_check(nc,search->second.getValues(&il1x_start));
       }
    }

    // Read number of frames.
    netcdf_check(nc,nframe = nc->ncid->getDim("number_of_images").getSize());
    check_error(nframe == 0,"Error: No frames in input file '%s'",filename.c_str());
    // Number of pixels in image is only relevant if there is an image.
    if (il1x_start < L1X_FOV) {
        netcdf_check(nc,npix = nc->ncid->getDim("samples_per_image").getSize());
    } else {
        npix = 0; // Irrelevant. If anyone uses this, it should be an error.
    }
    // Read metadata.
    NcGroup grp;
    // Group: image attributes.
    // Read reference attribute of image_time. This will be written in L1B.
    // For calibration L1A, you will have different files and possibly
    // inconsistent reference. Then, the L1B will raise a warning.
    // Operationally, L1B will be done with flight L1A, thus one file.
    // If you test L1B with calibration measurements, you can ignore the
    // warning, or you have files with same reference and then this will
    // be taken.

    { // Search variable attributes.
        map<string,NcVarAtt> mp;
        netcdf_check(nc,mp = nc->ncid->getGroup("image_attributes").getVar("image_time").getAtts());
        map<string,NcVarAtt>::iterator search;
        search = mp.find("reference");
        if (search == mp.end()) time_reference = "Unknown";
        else {
            netcdf_check(nc,search->second.getValues(time_reference));
        }
    }

    // Dimensions from science data (only for high L1X levels).
    if (
        il1x_start == L1X_UNBIN ||
        il1x_start == L1X_STRAY ||
        il1x_start == L1X_FOV ||
        il1x_start == L1X_RAD
    ) {
        netcdf_check(nc,grp = nc->ncid->getGroup("science_data"));
        if (
            il1x_start == L1X_UNBIN ||
            il1x_start == L1X_STRAY
        ) {
            // Spatial and spectral detector dimension.
            netcdf_check(nc,l1x_dim_detector_spat = grp.getDim("detector_spat").getSize());
            netcdf_check(nc,l1x_dim_detector_spec = grp.getDim("detector_spec").getSize());
        }
        // Number of spectra.
        if (
            il1x_start == L1X_FOV ||
            il1x_start == L1X_RAD
        ) {
            // Spatial detector dimension.
            netcdf_check(nc,l1x_dim_fov = grp.getDim("spatial_samples_per_image").getSize());
        }
    }

    // Close the file, because there may be an onslaught of L1A files.
    nc->close();

    return 0;

} // }}}

} // namespace tango
