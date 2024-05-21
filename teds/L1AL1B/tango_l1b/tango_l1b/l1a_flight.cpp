// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "array.h"
#include "settings_proc.h"
#include "l1a_file_metadata.h"
#include "l1a.h"
#include "l1a_flight.h"

namespace tango {

// Constructor.
L1A_flight::L1A_flight( // {{{
    Logger *creator
) : L1A_manager(creator)
{
} // }}}
// Virtual image counter.
int L1A_flight::count_images( // {{{
    size_t nfile, // Number of files.
    Settings_proc *set, // Processor settings.
    unique_ptr<L1A_file_metadata> *file_meta, // File-specific metadata, use only member 'nframe'.
    size_t &nl1a_total // Output: Number of L1A images in the entire run.
)
{

    // 1. Verify that there is only one file.
    // 2. Read the subset.
    // 3. If no subset is given, take the entire file.

    // 1. Verify that there is only one file.
    // For flight L1A, there should be only one L1A file in the settings.
    check_error(nfile != 1,"Error: Flight L1A works only for exactly one L1A file.");

    // 2. Read the subset size.
    check_error(set->image_start.size() != set->image_end.size(),"Error: Some selected image streaks either have no start or no end.");

    nl1a_total = set->image.size(); // Loose images.
    for (size_t istreak=0 ; istreak<set->image_start.size() ; istreak++) { // Loop over streaks.
        check_error(set->image_start[istreak] > set->image_end[istreak],"Error: Backward image streak from %zu to %zu\n",set->image_start[istreak],set->image_end[istreak]);
        nl1a_total += set->image_end[istreak] - set->image_start[istreak] + 1;
    }

    // 3. If no subset is given, take the entire file.
    if (nl1a_total == 0) nl1a_total = file_meta[0]->nframe;

    return 0;

} // }}}
// Virtual (subset-)organizer.
int L1A_flight::organize( // {{{
    size_t nfile, // Number of files.
    Settings_proc *set, // Processor settings.
    unique_ptr<L1A_file_metadata> *file_meta // File-specific metadata, use only member 'nframe' and set a pointer.
)
{

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

    // 1. Make sure each image has its shared pointer to the right NetCDF object.
    // There is only one file, so only one NetCDF object.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        instances_all[il1a]->file_meta = file_meta[0].get();
    }

    // 2. Assign the correct frames to each image.
    // 3. Verify that these frames exist.
    // It is already verified that each streak has a beginning and an end and that the end is
    // not before the beginning. The check for existence and overlap will
    // be done here.
    size_t nstreak = set->image_start.size();
    size_t nloose = set->image.size();
    if (nstreak + nloose == 0) {
        // 2. Assign the correct frames to each image.
        // No subset, so everything is included and no verification is needed.
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) { // nl1a_total is file_meta[0]->nframe.
            instances_all[il1a]->frame_offset = il1a;
        }
        // 3. Verify that these frames exist.
        // No action required.
    } else {
        // 2. Assign the correct frames to each image.
        // Define subset array, because it must be sorted and checked for overlap.
        vector<size_t> subset(nl1a_total);
        vector<size_t> subset_unsorted(nl1a_total);
        memcpy(subset_unsorted.data(),set->image.data(),nloose*sizeof(size_t));
        size_t *isubset_cur = &subset_unsorted[nloose];
        for (size_t istreak=0 ; istreak<nstreak ; istreak++) {
            for (size_t image=set->image_start[istreak] ; image<=set->image_end[istreak] ; image++) {
                *isubset_cur = image;
                isubset_cur++;
            }
        }
        Array<size_t>::sort(nl1a_total,subset_unsorted.data(),subset.data());
        // 3. Verify that these frames exist.
        // Negative indices need not be checked, because we have unsigned integers.
        check_error(subset[nl1a_total-1] >= file_meta[0]->nframe,"Error: Image index %zu not in file. Number of images in file is %zu.",subset[nl1a_total-1],file_meta[0]->nframe);
        // Verify for overlap.
        for (size_t il1a=1 ; il1a<nl1a_total ; il1a++) {
            check_error(subset[il1a] == subset[il1a-1],"Error: Image %zu is included more than once in the subset.",subset[il1a]);
        }
        // Put subset into L1A instances.
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) instances_all[il1a]->frame_offset = subset[il1a];
    }

    // All images cover only one frame, by definition. And in the L1X, they indexed in
    // order. Skipped images do not take fillvalue space to preserve frame indices. That
    // is not worth it.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        instances_all[il1a]->nframe_inp = 1;
        instances_all[il1a]->l1x_iframe = il1a;
    }
    // Warn if the settings desired a discard protocol. No discard protocol
    // is possible, because there is only one frame to start with.
    if (set->l1a_discard_frames_attempts != 0) {
        writelog(log_warning,"Warning: Frame discard protocol impossible for flight L1A. Ignoring discard settings.");
    }

    return 0;

} // }}}

} // namespace tango
