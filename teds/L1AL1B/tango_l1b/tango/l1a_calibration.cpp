// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings_proc.h"
#include "l1a_file_metadata.h"
#include "l1a.h"
#include "l1a_calibration.h"

// Constructor.
L1A_calibration::L1A_calibration( // {{{
    Logger *creator
) : L1A_manager(creator)
{
} // }}}
// Virtual image counter.
int L1A_calibration::count_images( // {{{
    size_t nfile, // Number of files.
    Settings_proc *set, // Processor settings.
    unique_ptr<L1A_file_metadata> *file_meta, // File-specific metadata, use only member 'nframe'.
    size_t &nl1a_total // Output: Number of L1A images in the entire run.
)
{

    // Each file is one image.
    nl1a_total = nfile;

    return 0;

} // }}}
// Virtual (subset-)organizer.
int L1A_calibration::organize( // {{{
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

    // Each file is one L1A image. Copy the relevant information from the files to the L1A instances.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        unique_ptr<L1A> &l1a = instances_all[il1a];
        l1a->file_meta = file_meta[il1a].get(); // LHS is a raw (or weak) pointer.
        // Perform the discard protocol.
        size_t discard_start = 0; // If all attempts fail, discard nothing.
        size_t discard_end = 0; // If all attempts fail, discard nothing.
        // Only allow frame discarding for pure L1A, not for L1X.
        if (file_meta[il1a]->il1x_start == L1X_L1A) {
            for (size_t iattempt=0 ; iattempt<set->l1a_discard_frames_attempts ; iattempt++) {
                size_t &at_start = set->l1a_discard_frames_start[iattempt];
                size_t &at_end = set->l1a_discard_frames_end[iattempt];
                if (at_start + at_end < file_meta[il1a]->nframe) {
                    discard_start = at_start;
                    discard_end = at_end;
                    break;
                }
            }
        }
        l1a->frame_offset = discard_start; // From the beginning, besides eventual discarded frames.
        l1a->nframe_inp = file_meta[il1a]->nframe - discard_start - discard_end; // The entire file, except for any discarded frames.
        l1a->l1x_iframe = 0; // There is only one frame.
    }

    return 0;

} // }}}

