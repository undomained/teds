// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "ckd.h"
#include "l1a.h"
#include "batch.h"
#include "processor.h"
#include "radcal.h"

// Settings functions.
Settings_radcal::Settings_radcal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "rad";
} // }}}
Settings_radcal::~Settings_radcal() {}
int Settings_radcal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // I think we need not do anything.
    // Perhaps later.
    return 0;

} // }}}

// Constructor for the Fovcal CKD structure.
Radcal::Radcal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("rad");
    set = make_unique<Settings_radcal>(this);
    Processor::set = set.get();
    own_skip = &ckd->rad_skip;
} // }}}
Radcal::~Radcal() {}

// Radiometric calibration software.
// This is the planned protocol.
// 1. Couple files to viewport.
// 2. Get radiometric calibration on spectrum level.
int Radcal::process_init( // {{{
)
{

    // Shape CKD.
    ckd->rad_spectra.resize(ckd->dim_fov*DIM_POL*ckd->dim_detector_spec,NC_FILL_DOUBLE);

    // Set viewport batches.
    ivp_batch.resize(ckd->dim_vp); // Maximum size.
    handle(batch_viewport(ivp_batch.data()));
    for (size_t ibatch=0 ; ibatch<nbatch ; ibatch++) {
        // Radiometric calibration only works with one image per viewport.
        check_error(batches[ibatch].nl1a != 1,"Error: Radiometric calibration needs exactly one image per viewport. Viewport %zu has %zu.\n",ivp_batch[ibatch],batches[ibatch].nl1a);
    }

    return 0;

} // }}}
int Radcal::process_batch(size_t ibatch, const Calibration_options& opt)
{

    // The viewport index is linked to the batch index.
    size_t &ivp = ivp_batch[ibatch];

    writelog(log_trace,"Calibrating radiance per spectrum.");

    double *rad_spectra_cur = &ckd->rad_spectra[ckd->dim_detector_spec]; // Running pointer.
    double *wave_spectra_cur = &ckd->wave_spectra[ckd->dim_detector_spec]; // Running pointer.

    L1A *l1a = l1a_instances[0];

    for (size_t ifov=0 ; ifov<ckd->fov_nfov_vp[ivp] ; ifov++) {

        Spectra specs;
        // Extract the spectrum.
        handle(l1a->extract(ifov, opt, specs));

        // Interpolate reference spectrum to extracted spectra.

        // So, we will use the wavelength domain of the extracted spectra.
        // There are two of them.
        // We will use the wavelengths from the spectrum-based wavelength calibration.
        // Actually, this does not even depend on the extracted spectra, just on the
        // spectra-based wavelength CKD.

        for (size_t ipol=0 ; ipol<DIM_POL ; ipol++) {
            vector<double> ref_spectrum(ckd->dim_detector_spec);
            linear_interpol(l1a->dim_refspec,ckd->dim_detector_spec,l1a->refspec_wavelength.data(),wave_spectra_cur,l1a->refspec_radiance.data(),ref_spectrum.data());
            // The intended radiance should be half the actual radiance, because this
            // is unpolarized light and then, S+ and S- should both be half the radiance.
            for (size_t ispec=0 ; ispec<ckd->dim_detector_spec ; ispec++) {
                rad_spectra_cur[ispec] = 0.5 * ref_spectrum[ispec] * l1a->exposure_time / specs.signal[ipol*ckd->dim_detector_spec+ispec];
                // Zeros, NaNs and infinites should never be used.
                if (!isnormal(rad_spectra_cur[ispec])) rad_spectra_cur[ispec] = NC_FILL_DOUBLE;
            }
            rad_spectra_cur += ckd->dim_detector_spec; // Move the pointer.
            wave_spectra_cur += ckd->dim_detector_spec; // Move the pointer.
        }
    }

    return 0;

} // }}}

