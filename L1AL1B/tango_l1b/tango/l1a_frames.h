// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1A_FRAMES_H
#define L1A_FRAMES_H

#include "header.h"
#include "l1a_manager.h"

class L1A_frames : public L1A_manager { // {{{

    public:
    L1A_frames(
        Logger *creator
    );
    private:
    int count_images(
        size_t nfile, // Number of files.
        Settings_proc *set, // Processor settings.
        unique_ptr<L1A_file_metadata> *file_meta, // File-specific metadata, use only member 'nframe'.
        size_t &nl1a_total // Output: Number of L1A images in the entire run.
    ) override;
    int organize(
        size_t nfile, // Number of files.
        Settings_proc *set, // Processor settings.
        unique_ptr<L1A_file_metadata> *file_meta // File-specific metadata, use only member 'nframe' and set a pointer.
    ) override;

    // Private administrative arrays.
    vector<size_t> il1a_start_file;
    vector<size_t> nl1a_file;
    vector<size_t> frame_offset_file;


}; // }}}

#endif
