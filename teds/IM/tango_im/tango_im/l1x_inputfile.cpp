#include "header.h"
#include "functions.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"
#include "l1x_inputfile.h"
#include "array.h"

namespace tango {

// Constructor.
L1X_inputfile::L1X_inputfile( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
// Destructor.
L1X_inputfile::~L1X_inputfile() {}

int L1X_inputfile::init( // {{{
    Settings_main *set,
    CKD *ckd
)
{
    nc = make_unique<NetCDF_object>(this);
    handle(nc->open(set->l1x_input,NcFile::read));

    multimap<string,NcGroupAtt> mp;
    netcdf_check(nc,mp = nc->ncid->getAtts());
    multimap<string,NcGroupAtt>::iterator search;
    search = mp.find("l1x_maturity_level");
    if (search == mp.end()) {
        il1x_start = L1X_TRUTH;
    } else {
        netcdf_check(nc,search->second.getValues(&il1x_start));
    }

    // Read number of frames.
    size_t nframe_total;
    nframe_total = nc->ncid->getDim("bins_along_track").getSize();
    check_error(nframe_total < 1,"Error: No frames in L1X input.");
    NcDim dimid_npix;
    size_t dimcheck;
    netcdf_check(nc,dimid_npix = nc->ncid->getDim("samples_per_image"));
    if (dimid_npix.isNull()) {
        check_error(il1x_start < L1X_FOV,"Error: Missing number of samples per image for L1X level that still involves detector images.");
    } else {
        netcdf_check(nc,dimcheck = dimid_npix.getSize());
        check_error(dimcheck != ckd->npix,"Error: Number of input pixels in L1X should correspond to the CKD for unbinned images in this program. Now, CKD desires %zu pixels and L1X file gives %zu.",ckd->npix,dimcheck);
    }

    // Read the subset size.
    check_error(set->image_start.size() != set->image_end.size(),"Error: Some selected image streaks either have no start or no end.");

    nframe = set->image.size(); // Loose images.
    for (size_t istreak=0 ; istreak<set->image_start.size() ; istreak++) { // Loop over streaks.
        check_error(set->image_start[istreak] > set->image_end[istreak],"Error: Backward image streak from %zu to %zu\n",set->image_start[istreak],set->image_end[istreak]);
        nframe += set->image_end[istreak] - set->image_start[istreak] + 1;
    }

    // 3. If no subset is given, take the entire file.
    bool full = nframe == 0; // Because we are ruining nframe.
    if (full) nframe = nframe_total;
    subset.resize(nframe);
    if (full) {
        // Assign the correct frames to each image.
        // No subset, so everything is included and no verification is needed.
        for (size_t iframe=0 ; iframe<nframe ; iframe++) {
            subset[iframe] = iframe;
        }
        // Verify that these frames exist.
        // No action required.
    } else {
        size_t nstreak = set->image_start.size();
        size_t nloose = set->image.size();
        // Assign the correct frames to each image.
        // Define subset array, because it must be sorted and checked for overlap.
        vector<size_t> subset_unsorted(nframe);
        memcpy(subset_unsorted.data(),set->image.data(),nloose*sizeof(size_t));
        size_t *isubset_cur = &subset_unsorted[nloose];
        for (size_t istreak=0 ; istreak<nstreak ; istreak++) {
            for (size_t image=set->image_start[istreak] ; image<=set->image_end[istreak] ; image++) {
                *isubset_cur = image;
                isubset_cur++;
            }
        }
        Array<size_t>::sort(nframe,subset_unsorted.data(),subset.data());
        // 3. Verify that these frames exist.
        // Negative indices need not be checked, because we have unsigned integers.
        check_error(subset[nframe-1] >= nframe_total,"Error: Image index %zu not in file. Number of images in file is %zu.",subset[nframe-1],nframe_total);
        // Verify for overlap.
        for (size_t iframe=1 ; iframe<nframe ; iframe++) {
            check_error(subset[iframe] == subset[iframe-1],"Error: Image %zu is included more than once in the subset.",subset[iframe]);
        }
    }

    NcGroup grp;
    // Dimensions from science data (only for high L1X levels).
    if (
        il1x_start == L1X_UNBIN ||
        il1x_start == L1X_STRAY ||
        il1x_start == L1X_FOV ||
        il1x_start == L1X_RAD ||
        il1x_start == L1X_TRUTH
    ) {
        netcdf_check(nc,grp = nc->ncid->getGroup("science_data"));
        if (
            il1x_start == L1X_UNBIN ||
            il1x_start == L1X_STRAY
        ) {
            // Spatial and spectral detector dimension.
            netcdf_check(nc,dimcheck = grp.getDim("detector_spat").getSize());
            check_error(dimcheck != ckd->dim_detector_spat,"Error: Inconistent spatial detector dimension between L1X (%zu) and CKD (%zu).",dimcheck,ckd->dim_detector_spat);
            netcdf_check(nc,dimcheck = grp.getDim("detector_spec").getSize());
            check_error(dimcheck != ckd->dim_detector_spec,"Error: Inconistent spectral detector dimension between L1X (%zu) and CKD (%zu).",dimcheck,ckd->dim_detector_spec);
        }
        // Number of spectra.
        if (
            il1x_start == L1X_FOV ||
            il1x_start == L1X_RAD
        ) {
            netcdf_check(nc,dimcheck = grp.getDim("spatial_samples_per_image").getSize());
            check_error(dimcheck != ckd->dim_fov,"Error: Inconistent number of spatial samples per image between L1X (%zu) and CKD (%zu).",dimcheck,ckd->dim_fov);
            netcdf_check(nc,dimcheck = grp.getDim("spectral_samples_per_spectrum").getSize());
            check_error(dimcheck != ckd->dim_detector_spec,"Error: Inconistent number of samples per spectrum between L1X (%zu) and CKD (%zu). Note that L1X spectra should be unbinned.",dimcheck,ckd->dim_detector_spec);
        }
        if (il1x_start == L1X_TRUTH) {

            // Read the dimensions of the input scene. Some have to comply
            // some CKD or hardcoded constants (viewport, Stokes parameters).
            // Others are arbitrary.
            size_t dimcheck;
            dim_spat_truth = nc->ncid->getDim("bins_across_track").getSize();
            // Read spectral dimension. It is a member variable, because it
            // will be needed during interpolation or ISRF convolution.
            dim_spec_truth = nc->ncid->getDim("bins_spectral").getSize();
            // Read wavelengths.
            wavelength.resize(ckd->dim_vp*dim_spec_truth);
            netcdf_check(nc,nc->ncid->getVar("wavelength").getVar(wavelength.data()));

        }
    }

    // GSE-data.
    // Group: gse_data (optional).
    netcdf_check(nc,grp = nc->ncid->getGroup("gse_data"));
    gse_exists = !grp.isNull();
    if (gse_exists) {
        { // Dimensions.
            multimap<string,NcDim> mp;
            netcdf_check(nc,mp = grp.getDims());
            multimap<string,NcDim>::iterator search;
            search = mp.find("wavelength");
            if (search == mp.end()) dim_refspec = 0;
            else {
                netcdf_check(nc,dim_refspec = search->second.getSize());
            }
            // Apply dimension.
            refspec_wavelength.resize(dim_refspec);
            refspec_radiance.resize(dim_refspec);
        }
        { // Variables.
            multimap<string,NcVar> mp;
            netcdf_check(nc,mp = grp.getVars());
            multimap<string,NcVar>::iterator search;
            search = mp.find("wavelength");
            if (search == mp.end()) {
                for (size_t iref=0 ; iref<dim_refspec ; iref++) refspec_wavelength[iref] = NC_FILL_DOUBLE; // Usually a zero-sized array.
            } else {
                netcdf_check(nc,search->second.getVar(refspec_wavelength.data()));
            }
            search = mp.find("signal");
            if (search == mp.end()) {
                for (size_t iref=0 ; iref<dim_refspec ; iref++) refspec_radiance[iref] = NC_FILL_DOUBLE; // Usually a zero-sized array.
            } else {
                netcdf_check(nc,search->second.getVar(refspec_radiance.data()));
            }
            search = mp.find("viewport");
            if (search == mp.end()) viewportinteger = NC_FILL_UBYTE;
            else {
                netcdf_check(nc,search->second.getVar(&viewportinteger));
            }
            
        }

        { // Attributes.
            multimap<string,NcGroupAtt> mp;
            netcdf_check(nc,mp = grp.getAtts());
            multimap<string,NcGroupAtt>::iterator search;
            search = mp.find("Illumination_level");
            if (search == mp.end()) illumination_level = NC_FILL_DOUBLE;
            else {
                netcdf_check(nc,search->second.getValues(&illumination_level));
            }
            search = mp.find("DoLP");
            if (search == mp.end()) refspec_dolp = NC_FILL_DOUBLE;
            else {
                netcdf_check(nc,search->second.getValues(&refspec_dolp));
            }
            search = mp.find("AoLP");
            if (search == mp.end()) refspec_aolp = NC_FILL_DOUBLE;
            else {
                netcdf_check(nc,search->second.getValues(&refspec_aolp));
            }
            search = mp.find("ACT_rotationAngle");
            if (search == mp.end()) act_angle = NC_FILL_DOUBLE;
            else {
                netcdf_check(nc,search->second.getValues(&act_angle));
            }
            search = mp.find("ALT_rotationAngle");
            if (search == mp.end()) alt_angle = NC_FILL_DOUBLE;
            else {
                netcdf_check(nc,search->second.getValues(&alt_angle));
            }
        }
    }

    return 0;

} // }}}

} // namespace tango
