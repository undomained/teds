// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "ckd.h"
#include "dimcal.h"

// Settings functions.
Settings_dimcal::Settings_dimcal( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "dim";
} // }}}
Settings_dimcal::~Settings_dimcal() {}
int Settings_dimcal::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(detector_spec); // Number of detector pixels in spectral dimension.
    recognize_setting(detector_spat); // Number of detector pixels in spatial dimension.
    recognize_setting(vp); // Number of viewports.

    return 0;

} // }}}

// Constructor for the Dimcal CKD structure.
Dimcal::Dimcal( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("dim");
    set = make_unique<Settings_dimcal>(this);
    Processor::set = set.get();
} // }}}
Dimcal::~Dimcal() {}

// Dimension calculation.
// This is the planned protocol.
// 1. Copy the dimensions from the settings.
// 2. Create empty pixel mask.
// 3. Create viewport mask.
int Dimcal::process_init( // {{{
)
{

    // Verify legal settings.
    check_error(set->detector_spec == 0,"Error: Missing or zero spectral detector dimension. Setting 'detector_spec'");
    check_error(set->detector_spat == 0,"Error: Missing or zero spatial detector dimension. Setting 'detector_spat'");
    check_error(set->vp == 0,"Error: Missing or zero viewports. Setting 'vp'");

    // 1. Copy the dimensions from the settings.
    ckd->dim_detector_spec = set->detector_spec;
    ckd->dim_detector_spat = set->detector_spat;
    ckd->dim_vp = set->vp;
    // Number of pixels in an image.
    ckd->npix = ckd->dim_detector_spec*ckd->dim_detector_spat;

    // 2. Create empty pixel mask.
    ckd->mask.resize(ckd->npix,false);
    ckd->vp_mask.resize(ckd->dim_vp,false);

    // 3. Create viewport mask.
    // Skip flagged viewports. The dimcal has to do this herself, because
    // the other processors do this before the process_init, so for dimcal,
    // the vp_mask is not yet created, because you do not yet know the
    // number of viewports.
    for (size_t ivp=0 ; ivp<set->vp_skip.size() ; ivp++) ckd->vp_mask[ivp] = set->vp_skip[ivp] == 1;

    return 0;

} // }}}

