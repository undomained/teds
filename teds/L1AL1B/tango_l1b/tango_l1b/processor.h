// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef PROCESSOR_H
#define PROCESSOR_H

#include "header.h"
#include "logger.h"
#include "l1a_manager.h" // For orbit_metadata.

class NetCDF_object;
class Settings_main;
class Settings_proc;
class CKD;
class L1A;
class Batch;

class Processor : public Logger { // {{{

    protected:
    // Protected constructor to avoid users from constructing a raw processor class.
    // The Processor should be considered as an abstract class.
    Processor(
        Logger *creator,
        CKD *ckd_arg
    );

    public:
    // Range of batches that each MPI process should process. The
    // starting batch index for a given process is mpi_ranges[my_rank]
    // and the end index is mpi_ranges[my_rank+1]-1. Example: for 7
    // images and 3 MPI processes we would have mpi_ranges { 0, 3, 5,
    // 7 } with the workloads 3,2,2. Each MPI process holds the full
    // vector, i.e. information about what the other processes are
    // doing. This is needed later in the synchronization step.
    std::vector<int> mpi_ranges {};

    // Virtual destructor, so if we want to use destructors in the end, we are
    // not facing nasty surprises in the future. This is because class Processor is
    // there to inherit from.
    virtual ~Processor();

    protected:
    // Pointer to the calibration key data.
    CKD *ckd;

    // Processor settings. This pointer will point to an instance of a child
    // of Settings_proc, made for the specific process.
    Settings_proc *set; // Processor settings structure (will point to derived member).

    // Orbit metadata.
    Orbit_metadata orb;

    // Array of all L1A instances, to be used in process_init, with just metadata.
    size_t nl1a_total;
    unique_ptr<L1A> *l1a_instances_all; // List of pointers to all L1A instances.
    // Batches, to be constructed in process_init, or kept as default.
    size_t nbatch;
    vector<Batch> batches;
    // List of L1A instances of the current batch, to be used during process_batch. Here, the image data is read and calibrated up to the desired level.
    size_t nl1a; // Number of L1A instances in current batch.
    L1A **l1a_instances; // List of L1A pointers in current batch.

    // Pointer to skipper.
    bool *own_skip = 0;

    // Member functions.
    public:
    // Main protocol of executing the processor.
    int execute(
        string &settings_file, // Settngs file name (from argument).
        Settings_main *set_main // Main settings, which are already read.
    );

    private:
    // Empty routines that need to be overwritten by the classes.
    // It is legal not to overwrite them, but then, nothing happens.
    virtual int process_init();
    virtual int process_batch(size_t ibatch, const Calibration_options& opt);
    virtual int process_finalize();
    virtual int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    );
    protected:
    // Standard batching algorithms.
    int batch_none();
    int batch_all();
    int batch_viewport(
        size_t *ivp_batch
    );
    int batch_one();

}; // }}}

#endif
