// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1A_MANAGER_H
#define L1A_MANAGER_H

#include "header.h"
#include "logger.h"

// Forward declaration.
class Calibration_options;
class Settings_proc;
class CKD;
class L1A_file_metadata;
class L1A;
class Batch;

// Struct for all metadata that has to do with the orbit. In principle,
// the L1B processor will just copy it to the output without looking at it.
// A CKD processor will probably ignore all this. And the L1A manager knows
// the content. They come from the files, but in principle, there should
// be only one file for L1B. However, if multiple files will ever be used,
// consistent metadata should just be copied and in the case of inconsistent
// metadata, a warning should be raised.
struct Orbit_metadata { // {{{

    // Global attributes to be copied from L1A without looking at them.
    string startDirection; // Ascending or descending.
    string endDirection; // Ascending or descending.
    long long orbit_number;
    string time_coverage_start;
    string time_coverage_stop;
    // Variable attribute to be copied from L1A without looking at it.
    string time_reference; // Reference on which image time is based.

}; // }}}

// L1A manager superclass.
class L1A_manager : public Logger { // {{{
    protected:
    // Protected constructor, so that no one constructs an undefined L1A manager.
    L1A_manager(
        Logger *creator
    );
    size_t nl1a_total; // Number of images, accessible for everyone.
    vector<unique_ptr<L1A>> instances_all; // All L1A-instances.

    public:
    virtual ~L1A_manager();
    // Getters.
    size_t getNtotal();
    unique_ptr<L1A> *getInstances();

    // The real contents.
    int init(
        Settings_proc *set,
        CKD *ckd,
        Orbit_metadata *orb
    );
    int l1x_init(
        Settings_proc *set,
        CKD *ckd
    );
    int construct_batch(
        CKD *ckd,
        Batch &bat,
        Calibration_options &opt,
        bool percentagelog_batch
    );
    int destruct_batch(
        Batch &bat
    );

    private:
    virtual int count_images(
        size_t nfile, // Number of files.
        Settings_proc *set, // Processor settings.
        unique_ptr<L1A_file_metadata> *file_meta, // File-specific metadata, use only member 'nframe'.
        size_t &nl1a_total // Output: Number of L1A images in the entire run.
    ) = 0;
    virtual int organize(
        size_t nfile, // Number of files.
        Settings_proc *set, // Processor settings.
        unique_ptr<L1A_file_metadata> *file_meta // File-specific metadata, use only member 'nframe' and set a pointer.
    ) = 0;

}; // }}}

#endif
