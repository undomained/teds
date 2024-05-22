// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "l1a_manager.h"
#include "l1a_flight.h"
#include "processor.h"

namespace tango {

// Constructor. This directly provides access to the CKD. Child-class
// constructors use the CKD, for the skip flag (if it is a skippable step).
// For the rest, child-class constructors will construct and link their own
// specific processor settings.
Processor::Processor( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Logger(creator)
{
    // The CKD should be available as soon as possible.
    ckd = ckd_arg;
} // }}}

Processor::~Processor() {}

// Executes the protocol for the execution of one actual run.
// One actual run is one process, with the settings of its section in
// the settings file.
int Processor::execute( // {{{
    string &settings_file,
    Settings_main *set_main
)
{

    // This is the protocol.
    // 1. Read the settings file.
    // 2. Figure out if the step is skipped.
    // 3. Construct the L1A manager.
    // 4. Perform initial actions of the processor before any L1A is imported.
    // 5. Loop over L1A batches.
    // 5a. Construct the batch.
    // 5b. Execute the processor per batch.
    // 6. Do whatever must be done at the end, after the last batch.
    // 7. Write the obtained CKD to output NetCDF file.
    // 8. Write detailed output if desired.

    // 1. Read the settings file.
    handle(set->init(settings_file));
    // Warn for non-existing output levels.
    if (set->output_level != NC_FILL_INT && set->output_level > set->output_level_max) {
        writelog(log_warning,"Warning: Requested output level (%d) higher than maximum supported by processor (%d).",set->output_level,set->output_level_max);
    }
    // Revert to default for omitted settings of which a default exists.
    handle(set->apply_main_settings(set_main));
    // Ignore any too high output level. If this is because of the
    // main settings, no warning is raised, but the result is the same.
    if (set->output_level > set->output_level_max) set->output_level = set->output_level_max;
    // Flag off viewports that are told to be skipped. If you have no CKD,
    // that is if you are L1C, this is no longer relevant. However, the L1C
    // processor will support flagging off viewports internally.
    if (ckd != NULL && ckd->lev > LEVEL_DIMCAL) { // Otherwise, the mask does not yet exist. The dimcal processor will do this after creation on its own.
        if (set->vp_skip.size() > 0) {
            for (size_t ivp=0 ; ivp<set->vp_skip.size() ; ivp++) {
                if (set->vp_skip[ivp] != 1) continue; // Ignore any 'false'. These are there if specific elements are set to 'true' in the settings.
                ckd->vp_mask[ivp] = true;
                // These settings only kill off viewports. You cannot revive any, because they will in general miss some CKD.
            }
        }
    }

    // 2. Figure out if the step is skipped.
    if (set->skip) {

        // Do not perform process at all.
        check_error(own_skip == 0,"Error: Attempt to skip a process that cannot be skipped.");
        writelog(log_info,"Skipping process.");
        set->output_level = 0; // Skip everything, also detailed output.
        *own_skip = true;

    } else {

        // Perform the process.
        writelog(log_info,"Starting process.");

        // Mark, if it is a skippable process, that it is not skipped.
        if (own_skip != 0) *own_skip = false;

        // Warn for consistency in optional calibration options.
        // Also, save current calibration options for later checks.
        // This is done at this position, because if the process itself
        // is skipped, nothing need to be checked.
        if (ckd != NULL) handle(ckd->check_opts(set->opt)); // L1C has nothing to do. It is only up to L1B.

        // 3. Construct the L1A manager.
        // Figure out how to batch the L1A images.

        // Prepare the images in a way, depending on what L1A style is chosen.
        // The settings is used to acquire information from, the L1A file list and
        // eventually some subset information.
        unique_ptr<L1A_manager> l1am;
        l1am = make_unique<L1A_flight>(this);

        handle(l1am->init(set,ckd,&orb)); // Constructor with error handling.

        // Smuggle the full L1A array from L1A manager to the processor.
        nl1a_total = l1am->getNtotal();
        l1a_instances_all = l1am->getInstances();


        // Default nbatch. If not set, revert to 'all' or 'none' based on nl1a_total.
        nbatch = NC_FILL_UINT64;

        // 4. Perform initial actions of the processor before any L1A is imported.
        // In this stage, also the batches can be configured. If no batch
        // configuration is executed, the program will revert to the default.

        handle(process_init()); // Before any batch is read. All L1A instances are available as l1a_instances_all, but only the metadata is read yet.
        // Revert to default batching if number of batches is still fillvalue.
        if (nbatch == NC_FILL_UINT64) {
            if (nl1a_total == 0) {
                handle(batch_none());
            } else {
                handle(batch_all());
            }
        }

        if (nbatch != 0) {
            writelog(log_info,"Number of batches: %zu",nbatch);
            writelog(log_info,"Total number of measurements: %zu",nl1a_total);
        }

        // Only write percentages if processor likes it. That is if there are
        // no percentage loggings in the processor itself.
        if (set->overall_percentages) {
            percentagelog_open("Processing images");
        }

        // 5. Loop over L1A batches.
        // In a parallel run the range of batches is given by
        // mpi_ranges. Otherwise the range is simply [0,nbatch).
        const size_t batch_begin {
            mpi_ranges.empty() ? 0 : static_cast<size_t>(mpi_ranges.at(my_rank))
        };
        const size_t batch_end {
            mpi_ranges.empty()
            ? nbatch
            : static_cast<size_t>(mpi_ranges.at(my_rank + 1))
        };
        for (size_t ibatch { batch_begin }; ibatch < batch_end; ibatch++) {
            // Write (update) percentage if desired.
            if (set->overall_percentages) {
                percentagelog_progress(ibatch, batch_end - batch_begin);
            }

            Batch &bat = batches[ibatch];

            // If ipix_end is not set, it will be interpreted to l1a->bin->npix,
            // where l1a is the L1A instance. This is done inside construct_batch
            // from the L1A manager.

            // 5a. Construct the batch.
            // This means that the L1A images (large data load) are read for the batch members.
            // And the detector calibration is performed.
            handle(l1am->construct_batch(ckd,bat,set->opt,!set->overall_percentages)); // The L1A images are owned by the L1A-manager.

            // The batch contains everything we need, but we extract the stuff we require for processing: the list of images and its size.
            nl1a = bat.nl1a;
            l1a_instances = bat.l1a.data();
            // 5b. Execute the processor per batch.
            handle(process_batch(ibatch, set->opt)); // Batch number may indicate the phase of the process.
            handle(l1am->destruct_batch(bat)); // Images are thrown away if not flagged for recycling. This is done to save memory.
        }

        if (set->overall_percentages) {
            percentagelog_close();
        }

        // 6. Do whatever must be done at the end, after the last batch.
        handle(process_finalize()); // After the last batch.
    }

    // 7. Write the obtained CKD to output NetCDF file.
    // This advances the CKD one step. This does not do anything for
    // the L1B or L1C processes, because they do not produce CKD.
    if (ckd != NULL) handle(ckd->writestep());

    // 8. Write detailed output if desired.
    // This is optional diagnostic output for CKD generation steps.
    // L1B and L1C should do this themselves, because this detailed output
    // is written to the CKD output file, which does only exist for the
    // CKD generation steps.
    if (my_rank == 0 && set->output_level >= 1) {
        // Create group in output.
        string groupname = format("detailed_output_%s",set->tag.c_str());
        NcGroup grp;
        netcdf_check(ckd->nc_ckd,grp = ckd->nc_ckd->ncid->addGroup(groupname.c_str()));
        // We pass the group and a reference to the NetCDF object for easier use.
        handle(write_detailed_output(ckd->nc_ckd.get(),grp));
        // Synchronize, so that an error in a later step does not destroy
        // all output.
        netcdf_check(ckd->nc_ckd,ckd->nc_ckd->ncid->sync());
    }

    return 0;

} // }}}

// Empty process routines. {{{
// These routines are virtual, but need not be overwritten. If you want that
// nothing is done, do not overwrite them.
int Processor::process_init(
)
{
    // Do nothing.
    return 0;
}
int Processor::process_batch(size_t ibatch, const Calibration_options& opt)
{
    // Do nothing.
    return 0;
}
int Processor::process_finalize(
)
{
    // Do nothing.
    return 0;
}
int Processor::write_detailed_output(
    NetCDF_object *nc,
    NcGroup &grp
)
{
    // Do nothing.
    return 0;
}
// }}}

// Standard batching routines. For more fancy batching, do it manually.
int Processor::batch_none( // {{{
)
{
    // No batch, no L1A. Nothing.
    nbatch = 0;
    batches.resize(nbatch);
    return 0;
} // }}}
int Processor::batch_all( // {{{
)
{
    // All L1A in one batch.
    nbatch = 1;
    batches.resize(nbatch);
    Batch &bat = batches[0];
    bat.setSize(nl1a_total);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        bat.l1a[il1a] = l1a_instances_all[il1a].get(); // From list of unique_ptrs to list of pointers.
    }
    return 0;
} // }}}
int Processor::batch_viewport( // {{{
    size_t *ivp_batch
)
{
    // Rather than knowing how many viewports there are, just look at the bits of the viewport array.
    // Get maximum viewport number for this.
    vector<uint8_t> viewportarray(nl1a_total);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) viewportarray[il1a] = l1a_instances_all[il1a]->viewport;
    nbatch = 0;
    for (size_t ivp=0 ; ivp<ckd->dim_vp ; ivp++) if (!ckd->vp_mask[ivp]) nbatch++;
    batches.resize(nbatch);
    size_t ivp = 0;
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        while (ckd->vp_mask[ivp]) {
            // Go to next viewport.
            for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
                viewportarray[il1a] /= 2; // Go to next view port.
            }
            ivp++;
        }
        if (ivp_batch != NULL) ivp_batch[ibatch] = ivp;
        Batch &bat = batches[ibatch];
        // First count.
        size_t nl1a_batch = 0;
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            if ((viewportarray[il1a] & 1) == 1) nl1a_batch++;
        }
        bat.setSize(nl1a_batch);
        // Now, do the processing.
        size_t cnt = 0; // Element counter.
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
            if ((viewportarray[il1a] & 1) == 1) {
                bat.l1a[cnt] = l1a_instances_all[il1a].get(); // Copy this pointer.
                bat.remove[cnt] = viewportarray[il1a] == 1; // Remove if this is the last viewport.
                cnt++;
            }
            viewportarray[il1a] /= 2; // Go to next view port.
        }
        ivp++;
    }
    return 0;
} // }}}
int Processor::batch_one( // {{{
)
{
    // One file at a time.
    nbatch = nl1a_total;
    batches.resize(nbatch);
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        Batch &bat = batches[ibatch];
        bat.setSize(1);
        bat.l1a[0] = l1a_instances_all[ibatch].get();
    }
    return 0;
} // }}}

} // namespace tango
