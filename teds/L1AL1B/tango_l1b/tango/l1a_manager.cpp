// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "logger.h"
#include "settings_proc.h"
#include "ckd.h"
#include "l1a_file_metadata.h"
#include "binningtable.h"
#include "l1x.h"
#include "l1a.h"
#include "batch.h"
#include "l1a_manager.h"

// Constructor.
L1A_manager::L1A_manager( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
L1A_manager::~L1A_manager() {} // Destructor (virtual).

// Getters.
size_t L1A_manager::getNtotal( // {{{
)
{
    return nl1a_total;
} // }}}
unique_ptr<L1A> *L1A_manager::getInstances( // {{{
)
{
    return instances_all.data();
}
// }}}

// This plays the role of the constructor, but it return an integer for error
// handling. This is a non-virtual routine that calls virtual routines if there
// are different implementation per type of L1A manager.
int L1A_manager::init( // {{{
    Settings_proc *set,
    CKD *ckd,
    Orbit_metadata *orb
)
{
    size_t nfile = set->l1a_files.size();
    // Make an L1A-file structure for each of the files.
    vector<unique_ptr<L1A_file_metadata>> file_meta(nfile);
    // Read all the metadata that is file-specific (low data load, hopefully).
    if (nfile != 0) percentagelog_open("Reading L1A metadata");
    for (size_t ifile=0 ; ifile<nfile ; ifile++) {
        percentagelog_progress(ifile,nfile);
        file_meta[ifile] = make_unique<L1A_file_metadata>(this);
        file_meta[ifile]->file_id = ifile;
        file_meta[ifile]->filename = set->l1a_files[ifile];
        file_meta[ifile]->l1a_navigation = set->l1a_navigation;
        handle(file_meta[ifile]->read_metadata());
    }
    if (nfile != 0) percentagelog_close();

    // Interpret orbit metadata. {{{
    if (nfile == 0) {
        orb->startDirection = "Unknown";
        orb->endDirection = "Unknown";
        orb->orbit_number = NC_FILL_INT64;
        orb->time_coverage_start = "Unknown";
        orb->time_coverage_stop = "Unknown";
        orb->time_reference = "Unknown";
    } else {
        orb->startDirection = file_meta[0]->startDirection;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->startDirection.compare(file_meta[ifile]->startDirection) != 0) {
                orb->startDirection = "Inconsistent";
                break;
            }
        }
        orb->endDirection = file_meta[0]->endDirection;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->endDirection.compare(file_meta[ifile]->endDirection) != 0) {
                orb->endDirection = "Inconsistent";
                break;
            }
        }
        orb->orbit_number = file_meta[0]->orbit_number;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->orbit_number != file_meta[ifile]->orbit_number) {
                orb->orbit_number = NC_FILL_INT64; // No distinction between unknown and inconsistent, unfortunately. There is only one fill value.
                break;
            }
        }
        orb->time_coverage_start = file_meta[0]->time_coverage_start;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->time_coverage_start.compare(file_meta[ifile]->time_coverage_start) != 0) {
                orb->time_coverage_start = "Inconsistent";
                break;
            }
        }
        orb->time_coverage_stop = file_meta[0]->time_coverage_stop;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->time_coverage_stop.compare(file_meta[ifile]->time_coverage_stop) != 0) {
                orb->time_coverage_stop = "Inconsistent";
                break;
            }
        }
        orb->time_reference = file_meta[0]->time_reference;
        for (size_t ifile=1 ; ifile<nfile ; ifile++) {
            if (orb->time_reference.compare(file_meta[ifile]->time_reference) != 0) {
                orb->time_reference = "Inconsistent";
                break;
            }
        }
    } // }}}

    // Interpret and check memory-effect settings.
    // A pity that we are writing into the settings, but otherwise, we
    // need a separate routine.
    set->l1a_discard_frames_attempts = set->l1a_discard_frames_start.size();
    size_t dimcheck = set->l1a_discard_frames_end.size();
    check_error(dimcheck != set->l1a_discard_frames_attempts,"Error: Number of attempts to mitigate memory effect must be the same for the beginning (%zu) and the end (%zu).",set->l1a_discard_frames_attempts,dimcheck);

    // Derived routine: Count L1A images.
    handle(count_images(nfile,set,file_meta.data(),nl1a_total));

    // Allocate the L1A-structures.
    instances_all.resize(nl1a_total);
    // Apparently, you cannot use make_unique inside resize.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) instances_all[il1a] = make_unique<L1A>(this);

    // Perform some organization. The tasks of the derived class here is the following.
    // 1. Make sure each image has its shared pointer to the right file object.
    // 2. Assign the correct frames to each image.
    // 3. Verify that these frames exist.
    // 4. Set frame index for L1X output.
    // Step four is also executed if there is no L1x output, because it is not so much work.
    // The following fields should be assigned.
    // - instances_all[...]->file_meta (L1A_file_metadata *)
    // - instances_all[...]->frame_offset (size_t)
    // - instances_all[...]->nframe_inp (size_t)
    // - instances_all[...]->l1x_iframe (size_t)
    // It is highly desirable not to use any netcdf actions. Use file_meta only to set a pointer
    // to the right file and to peak into the member 'nframe'.
    handle(organize(nfile,set,file_meta.data()));

    // Read the metadata of all the L1A-instances.
    if (nl1a_total != 0) writelog(log_trace,"Organizing L1A metadata.");
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        unique_ptr<L1A> &l1a = instances_all[il1a];
        // Copy request for navigation data.
        l1a->l1a_navigation = set->l1a_navigation;
        // Create array of L1X pointers for L1X output.
        l1a->l1x.resize(nl1x); // This is an array of pointers. Nothing serious is allocated yet.
        // Read the metadata.
        handle(l1a->read_metadata());
    }

    // Investigate the binning tables.
    bool bin_file_needed = false; // Only needed if any non-trivial binning table is used.
    vector<uint8_t> binningtables(nl1a_total); // Maximum size.
    size_t nuniq = 0;
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        uint8_t &bin_cur = instances_all[il1a]->binning_table_id;
        bool unique = true;
        for (size_t iuniq=0 ; iuniq<nuniq ; iuniq++) {
            if (binningtables[iuniq] == bin_cur) {
                unique = false;
                break;
            }
        }
        if (unique) {
            binningtables[nuniq] = bin_cur;
            nuniq++;
            if (bin_cur != 0) bin_file_needed = true;
        }
    }

    // Clipping the binningtables array is not needed, because we have nuniq all the time.
    // Read the binning tables.
    unique_ptr<Binningtable_file> bin_file;
    if (bin_file_needed) {
        writelog(log_trace,"Constructing binning table.");
        bin_file = make_unique<Binningtable_file>(this);
        handle(bin_file->read(set->binningtable_filename,ckd));
    }
    // For L1X, keep track of largest spectrum size (only when FOV CKD already exists).
    size_t largest_spectrum_size = 0;
    for (size_t iuniq=0 ; iuniq<nuniq ; iuniq++) {
        uint8_t bin_cur = binningtables[iuniq];
        shared_ptr<Binningtable> table_cur = make_shared<Binningtable>(this);
        handle(table_cur->read(bin_file.get(),bin_cur,ckd)); // If this is identifier 0 and there is no file, the null pointer bin_file.get() suffices.
        // Distribute.
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            if (instances_all[il1a]->binning_table_id == bin_cur) {
                // Copy the shared pointer. This makes the object survive outside the loop.
                instances_all[il1a]->bin = table_cur;
            }
        }
        // Measure size of extracted spectra.
        if (ckd->lev > LEVEL_FOVCAL) {
            for (size_t ifov=0 ; ifov<ckd->dim_fov ; ifov++) {
                size_t &spectrumsize_cur = table_cur->binned_ckd->fov_dims_spec[ifov];
                if (spectrumsize_cur > largest_spectrum_size) largest_spectrum_size = spectrumsize_cur;
            }
        }
    }
    // Verify if the number of pixels is sufficient for the binning tables.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        L1A *l1a = instances_all[il1a].get();
        check_error(l1a->file_meta->il1x_start < L1X_FOV && l1a->file_meta->npix < l1a->bin->npix,"Error: File '%s' only has %zu pixels per image, but it includes binning table %d, which requires %zu pixels.",l1a->filename.c_str(),l1a->file_meta->npix,l1a->binning_table_id,l1a->bin->npix);
    }

    // Perform L1X administration.
    // The user sets L1X output by defining the output files per L1X
    // output step.
    vector<vector<string> *> l1x_outputfiles = {
        &set->l1x_outputfiles_raw,
        &set->l1x_outputfiles_dark,
        &set->l1x_outputfiles_noise,
        &set->l1x_outputfiles_nonlin,
        &set->l1x_outputfiles_prnu,
        &set->l1x_outputfiles_unbin,
        &set->l1x_outputfiles_stray,
        &set->l1x_outputfiles_rebin,
        &set->l1x_outputfiles_fov,
        &set->l1x_outputfiles_rad
    };

    // First count the images per file using the file_id property
    // from the L1A instances, which is metadata. This is required if
    // any L1X instance is written. And the same counts are relevant
    // for all L1X. Pity for the several nanoseconds calculating this
    // if no L1X output is selected.
    vector<size_t> l1x_nframe_file(nfile,0);
    // Each L1A instances increments its file identifier's index.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) l1x_nframe_file[instances_all[il1a]->file_id]++;
    for (size_t il1x=0 ; il1x<nl1x ; il1x++) {
        // Make an instance per file.
        vector<shared_ptr<L1X>> l1x(nfile);
        // Apparently, you cannot use make_shared in vector constructor.
        for (size_t ifile=0 ; ifile<nfile ; ifile++) l1x[ifile] = make_shared<L1X>(this);
        size_t sz = l1x_outputfiles[il1x]->size();
        if (sz == 0) {
            // No L1X chosen.
            // Important to have this clause here (or use another construction). If the
            // sz is equal to zero, the zero should be recognized. Important for the L1A-less
            // calibration step 'dim'.
            // Do nothing.
        } else if (sz == nfile) {
            // Yes. L1X chosen. Make the L1X instances mature.
            for (size_t ifile=0 ; ifile<nfile ; ifile++) {
                handle(l1x[ifile]->init(
                    (l1x_t) il1x,
                    l1x_nframe_file[ifile],
                    (*l1x_outputfiles[il1x])[ifile],
                    file_meta[ifile].get(),
                    ckd,
                    largest_spectrum_size
                ));
            }
        } else {
            // Wrong input.
            raise_error("Error: Number of L1X files for step %s is %zu, which is neither zero nor %zu, the number of L1A files.",l1x_steps[il1x].c_str(),sz,nfile);
        }
        // Distribute the L1X instances over the L1A instances.
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            L1A *l1a = instances_all[il1a].get();
            l1a->l1x[il1x] = l1x[l1a->file_id]; // Copy of shared pointer.
            handle(l1a->l1x[il1x]->write_metadata(l1a));
            // It is possible (likely for flight L1A) that multiple L1A instances get the
            // same L1X pointer.
        }
    }

    return 0;

} // }}}

// Construct and cleans up the batch.
int L1A_manager::construct_batch( // {{{
    CKD *ckd,
    Batch &bat,
    Calibration_options &opt,
    bool percentagelog_batch
)
{
    if (percentagelog_batch) percentagelog_open("Reading batch");
    for (size_t il1a=0 ; il1a<bat.nl1a ; il1a++) {
        if (percentagelog_batch) percentagelog_progress(il1a,bat.nl1a);

        L1A *&l1a = bat.l1a[il1a];
        // If the image is already there, do nothing.
        size_t bat_ipix_end = bat.ipix_end == NC_FILL_UINT64?l1a->bin->npix:bat.ipix_end;
        if (bat.ipix_start < l1a->ipix_start || bat_ipix_end > l1a->ipix_end) {

            // Read the high-data-load contents.
            // A fraction of zero is interpreted as the total number of pixels.
            handle(l1a->read_image(bat.ipix_start,bat_ipix_end));
            // Execute detector calibration.
            handle(l1a->calibrate_detector(ckd,opt));

        } else {

            writelog(log_verbose,"Recycling image from '%s'.",l1a->filename.c_str());

        }

    }
    if (percentagelog_batch) percentagelog_close();

    return 0;

} // }}}
int L1A_manager::destruct_batch( // {{{
    Batch &bat
)
{
    // Remove image if needed.
    for (size_t il1a=0 ; il1a<bat.nl1a ; il1a++) {
        if (bat.remove[il1a]) {
            handle(bat.l1a[il1a]->remove_image());
        }
    }
    return 0;

} // }}}

