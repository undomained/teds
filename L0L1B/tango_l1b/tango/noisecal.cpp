// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "matrix.h"
#include "lininv.h"
#include "netcdf_object.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "noisecal.h"

// Settings functions.
Settings_noisecal::Settings_noisecal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "noise";
    output_level_max = 1;
    l1a_type = L1A_FRAMES;
    opt.dark_current = false;
} // }}}
Settings_noisecal::~Settings_noisecal() {}
int Settings_noisecal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(nfrac); // Number of fractions to cut the L1A images.
    recognize_setting(mask_g_min); // Minimum signal-dependent noise term.
    recognize_setting(mask_g_max); // Maximum signal-dependent noise term.
    recognize_setting(mask_n_min); // Minimum signal-independent noise term.
    recognize_setting(mask_n_max); // Maximum signal-independent noise term.

    return 0;

} // }}}

// Constructor for the Noisecal CKD structure.
Noisecal::Noisecal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("noise");
    set = make_unique<Settings_noisecal>(this);
    Processor::set = set.get();
    own_skip = &ckd->noise_skip;
} // }}}
Noisecal::~Noisecal() {}

int Noisecal::process_init( // {{{
)
{

    // Shape target CKD.
    ckd->noise_g.resize(ckd->npix);
    ckd->noise_n.resize(ckd->npix);

    // Verify that all L1A have non-pre-coadded frames.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        L1A *l1a = l1a_instances_all[il1a].get();
        check_error(l1a->nr_coadditions != 1,"Error: Noise calibration only works with frames with co-addition factor one. File '%s' shows coaddition factor %d found at frame %zu.",l1a->filename.c_str(),l1a->nr_coadditions,l1a->frame_offset);
    }

    // Make a batch per fraction, with all L1A inside, but with a fraction
    // of the pixels.
    nbatch = set->nfrac;
    batches.resize(nbatch);
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        Batch &bat = batches[ibatch];
        bat.setSize(nl1a_total);
        for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) bat.l1a[il1a] = l1a_instances_all[il1a].get(); // All L1A instances.
        bat.ipix_start = ckd->npix * ibatch / nbatch;
        bat.ipix_end = ckd->npix * (ibatch+1) / nbatch;
    }

    // Get number of L1A files. Use the fact that the files are in order.
    // So the last L1A instance is from the last file.
    nfile = l1a_instances_all[nl1a_total-1]->file_id + 1;
    nl1a_file.resize(nfile,0);
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        nl1a_file[l1a_instances_all[il1a]->file_id]++;
    }
    // Loop should be over files, but L1A index is needed to clearly formulate the error message.
    // Unfortunately, this means there will be some double checks.
    // I hope this costs not too much time.
    for (size_t il1a=0 ; il1a<nl1a_total ; il1a++) {
        check_error(nl1a_file[l1a_instances_all[il1a]->file_id] < 2,"Error: File '%s' has less than two frames, so cannot be used to estimate noise.",l1a_instances_all[il1a]->filename.c_str());
    }

    // Shape optional output.
    if (set->output_level >= 1) {
        det_signal.resize(ckd->npix*nfile);
        det_noise.resize(ckd->npix*nfile);
    }

    return 0;

} // }}}

// Noise calibration protocol.
int Noisecal::process_batch( // {{{
    size_t ibatch
)
{

    // Get pixel range from batch.
    Batch &bat = batches[ibatch];
    // Now, the pixel range is bat.ipix_start to bat.ipix_end.

    // The formula that we fit (from ATBD) is:
    // Noise = sqrt(g*signal + n^2)
    // That is:
    // Noise^2 = g*signal + n^2
    // and Noise^2 = 1/(N-1) * sum((x-avg)^2)

    // We will fit noise^2 and the state parameters are n^2 and g.
    // The noise^2 is the measurement.
    // The Jacobian for g is the signal, for n^2 is one.

    // Running CKD pointers.
    double *noise_g_cur = &ckd->noise_g[bat.ipix_start];
    double *noise_n_cur = &ckd->noise_n[bat.ipix_start];

    // First, default to standard deviation, because we do not expect
    // outliers. Also, we expect limited data sets.
    percentagelog_open("Fitting noise parameters per pixel");
    for (size_t ipix=bat.ipix_start ; ipix<bat.ipix_end ; ipix++) {

        percentagelog_progress(ipix,ckd->npix);

        vector<double> jac(2*nfile,1.0); // The first half will be the average signal. The second half will stay one.
        vector<double> meas(nfile); // These will be the variances.
        vector<double> res(2);

        // Not necessary, but desirable to avoid naming confusion.
        vector<double> agg(nfile,0.0); // Aggregate.
        vector<double> aggsq(nfile,0.0); // Aggregate squares.
        for (size_t il1a=0 ; il1a<nl1a ; il1a++) {
            L1A *l1a = l1a_instances[il1a];
            size_t &ifile = l1a->file_id;
            double &dat = l1a->image[ipix];
            agg[ifile] += dat;
            aggsq[ifile] += pow(dat,2.0);
        }
        // Calculate average and variance.
        for (size_t ifile=0 ; ifile<nfile ; ifile++) {
            double &avg = jac[ifile];
            double &var = meas[ifile];
            avg = agg[ifile] / nl1a_file[ifile];
            var = (aggsq[ifile]/nl1a_file[ifile] - pow(avg,2.0)) * nl1a_file[ifile] / (nl1a_file[ifile]-1.0);
        }
        if (set->output_level >= 1) {
            // Write raw results into detailed output.
            memcpy(&det_signal[ipix*nfile],jac.data(),nfile*sizeof(double)); // Copy first nfile elements of jac. Those are the signals.
            // Noise needs another square root. All variances should be positive, because they are actual variances.
            // If they are zero or slightly negative through rounding errors, set noise to zero, possibly for a dead pixel or so.
            for (size_t ifile=0 ; ifile<nfile ; ifile++) det_noise[ipix*nfile+ifile] = meas[ifile]<=0.0?0.0:sqrt(meas[ifile]); // The variances are in meas.
        }
        check_error(linear_invert(2,nfile,jac.data(),OPT_NONE,NULL,NULL,meas.data(),res.data()),"Error: Noise parameter fit failed.");
        *noise_g_cur = res[0];
        // Guard for negative N squared.
        if (res[1] < 0.0) {
            *noise_n_cur = NC_FILL_DOUBLE;
            ckd->mask[ipix] = true;
        } else *noise_n_cur = sqrt(res[1]);

        // Apply pixel mask criteria.
        if (set->mask_g_min != NC_FILL_DOUBLE && *noise_g_cur < set->mask_g_min) ckd->mask[ipix] = true;
        if (set->mask_g_max != NC_FILL_DOUBLE && *noise_g_cur > set->mask_g_max) ckd->mask[ipix] = true;
        if (set->mask_n_min != NC_FILL_DOUBLE && *noise_n_cur < set->mask_n_min) ckd->mask[ipix] = true;
        if (set->mask_n_max != NC_FILL_DOUBLE && *noise_n_cur > set->mask_n_max) ckd->mask[ipix] = true;

        // Progress pointers.
        noise_g_cur++;
        noise_n_cur++;

    }
    int percentage_end = 100*bat.ipix_end/ckd->npix;
    percentagelog_close(percentage_end);

    return 0;

} // }}}

int Noisecal::write_detailed_output( // {{{
    NetCDF_object *nc,
    NcGroup &grp
)
{
    NcDim dimid_nmeas;
    netcdf_check(nc,dimid_nmeas = grp.addDim("nmeas",nfile));
    // Type-cast the size_t into int.
    vector<int> int_nl1a_file(nfile);
    for (size_t ifile=0 ; ifile<nfile ; ifile++) int_nl1a_file[ifile] = (int) nl1a_file[ifile];
    netcdf_check(nc,grp.addVar("nframe",ncInt,dimid_nmeas).putVar(int_nl1a_file.data()));
    vector<NcDim> dims = {ckd->dimid_detector_spat,ckd->dimid_detector_spec,dimid_nmeas};
    netcdf_check(nc,grp.addVar("signal",ncDouble,dims).putVar(det_signal.data()));
    netcdf_check(nc,grp.addVar("noise",ncDouble,dims).putVar(det_noise.data()));

    return 0;

} // }}}

