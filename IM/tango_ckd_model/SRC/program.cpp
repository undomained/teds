#include "functions.h"
#include "header.h"
#include "git.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "settings_l1b.h"
#include "settings_isrf.h"
#include "settings_noise.h"
#include "settings_geo.h"
#include "planet.h"
#include "ckd.h"
#include "l1x_inputfile.h"
#include "frame.h"
#include "l1x.h"
#include "binningtable.h"
#include "utc.h"
#include "program.h"
#include <random>

Program::Program( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}

int Program::execute( // {{{
    string settings_file,
    bool foldsettings
)
{

    writelog(log_info,"%s %s (%s)", "CKD instrument model", git_tag.c_str(), git_date.c_str());
    writelog(log_debug,"Commit: %s",git_commit.c_str());

    // Log the settings into the log file.
    // This scope makes sure to clean any I/O-related stuff from the stack.
    // It would remain there until the end of the program.
    {
        struct stat settingsfile_info;
        check_error(stat(settings_file.c_str(),&settingsfile_info),"Error acquiring file information from '%s'.",settings_file.c_str());
        string modtime = ctime(&settingsfile_info.st_mtime);
        size_t found = modtime.find("\n");
        modtime = modtime.substr(0,found);
        writelog(log_debug,"Settings file: '%s' (Last modified %s).%s",settings_file.c_str(),modtime.c_str(),foldsettings?" {{{":"");
        ifstream setstream(settings_file.c_str());
        check_error(setstream.fail(),"Error opening settings file: '%s'.",settings_file.c_str());
        string line;
        while (getline(setstream,line)) {
            writelog(log_debug,line.c_str());
        }
        writelog(log_debug,"End of settings file.%s",foldsettings?" }}}":"");
    }

    // Read settings.
    Settings_main set_main(this);
    handle(set_main.init(settings_file));

    // Construct and read the CKD.
    unique_ptr<CKD> ckd = make_unique<CKD>(this);
    handle(ckd->read(set_main,LEVEL_L1B,false));
    if (!set_main.dark_apply) {
        ckd->dark_skip = true;
    }
    if (!set_main.nonlin_apply) {
        ckd->nonlin_skip = true;
    }
    if (!set_main.prnu_apply) {
        ckd->prnu_skip = true;
    }
    if (!set_main.stray_apply) {
        ckd->stray_skip = true;
    }
    if (!set_main.rad_apply) {
        ckd->rad_skip = true;
    }

    size_t il1x_start;
    L1X_inputfile l1x_inputfile(this);

    handle(l1x_inputfile.init(&set_main,ckd.get()));
    il1x_start = l1x_inputfile.il1x_start;
    size_t nframe = l1x_inputfile.nframe; // Number of frames determined by L1X input.

    Settings_noise set_ns(this);
    set_ns.init(settings_file); // Always from input.
    Settings_l1b set_l1b(this);
    Settings_geo set_geo(this);
    unique_ptr<Planet> earth;
    unique_ptr<UTC> utc;
    if (il1x_start == L1X_TRUTH) {
        handle(set_l1b.init(settings_file));
    } else {
        set_l1b.execute = false;
    }

    // Construct all L1X.
    vector<unique_ptr<L1X>> l1x(il1x_start);
    for (size_t il1x=0 ; il1x<il1x_start ; il1x++) l1x[il1x] = make_unique<L1X>(this);
    vector<string *> l1x_outputfiles = {
        &set_main.l1x_outputfile_raw,
        &set_main.l1x_outputfile_dark,
        &set_main.l1x_outputfile_noise,
        &set_main.l1x_outputfile_nonlin,
        &set_main.l1x_outputfile_prnu,
        &set_main.l1x_outputfile_unbin,
        &set_main.l1x_outputfile_stray,
        &set_main.l1x_outputfile_rebin,
        &set_main.l1x_outputfile_fov,
        &set_main.l1x_outputfile_rad
    }; // Only using the first il1x_start elements.

    // Never write L1X from il1x_start or more towards the L1B.
    int il1x_end = nl1x;
    for (size_t il1x=0 ; il1x<il1x_start ; il1x++) {
        if (l1x_outputfiles[il1x]->compare("") != 0) {
            handle(l1x[il1x]->init(
                static_cast<l1x_t>(il1x), // L1X level identifier.
                nframe, // Number of frames.
                *l1x_outputfiles[il1x], // L1X file name.
                ckd->npix, // Full-detector number of pixels,
                0, // No binning table.
                ckd.get(), // Calibration key data (for detector and spectrum dimensions).
                &l1x_inputfile // For GSE-data to write.
            ));
            if (il1x_end == nl1x) il1x_end = il1x;
        }
    }
    unique_ptr<L1X> l1a = make_unique<L1X>(this);

    size_t npix_binned;
    unique_ptr<Binningtable> bin = make_unique<Binningtable>(this);
    if (set_main.binning_table_id == 0) {
        npix_binned = ckd->npix;
    } else {
        unique_ptr<Binningtable_file> bin_file = make_unique<Binningtable_file>(this);
        handle(bin_file->read(set_main.binningtable_filename,ckd->dim_detector_spat,ckd->dim_detector_spec));
        handle(bin->read(bin_file.get(),static_cast<uint8_t>(set_main.binning_table_id)));
        npix_binned = bin->npix;
    }

    if (set_main.l1a_outputfile.compare("") != 0) {
        handle(l1a->init(
            L1X_L1A,
            nframe,
            set_main.l1a_outputfile,
            npix_binned,
            set_main.binning_table_id,
            ckd.get(),
            &l1x_inputfile
        ));
        il1x_end = -1; // L1X_L1A.
    }

    for (size_t iframe=0 ; iframe<nframe ; iframe++) {
        std::cout << "PROCESSING FRAME " << iframe << std::endl;

        unique_ptr<Frame> frm = make_unique<Frame>(this);

        // Read L1X input.
        handle(frm->read_l1x(
            iframe,
            ckd.get(),
            &l1x_inputfile
        ));
        frm->exposure_time = set_main.exposure_time;
        frm->nr_coadditions = set_main.nr_coadditions;

        if (il1x_start == L1X_TRUTH) {

            // Add bonus options to truth.

            // Enhance input spectral sampling and interpolate to
            // CKD spatial samples. The spectral samples are not yet
            // CKD, because an ISRF convolution still be executed.
            // The spatial resampling is necessary, because we will construct
            // a reference truth L1B file from the finest possible spectral
            // samplig in combination with the CKD spatial sampling.
            // Enhancing the spectral sampling may improve the ISRF convolution
            // as well as the truth L1B gauss convolution.
            // Note that the interpolation will be done in the input units,
            // so before unit conversion.
            handle(frm->resample(ckd.get(),set_main.input_enhance));
            if (ckd->fov_nfov_vp.front() != frm->dim_spat_truth) {
                std::cerr << "error: the input file contains "
                          << frm->dim_spat_truth
                          << " spectra whereas the CKD defines "
                          << ckd->fov_nfov_vp.front() << " ACT angles\n";
                return 1;
            }

            // Convert photons to Joules. The rest of the units are
            // normal. There are no centimeters rubbish like that.
            if (set_main.photons) {
                handle(frm->convert_units(ckd.get()));
            }
            // Read ISRF settings including the flag whether or not there
            // is an ISRF.
            Settings_isrf set_isrf(this);
            handle(set_isrf.init(settings_file));
            if (set_isrf.execute) {
                frm->isrf_interpolate(ckd.get(), &set_isrf);
            } else {
                // Interpolate on CKD spectral resolution
                handle(frm->interpolate_truth(ckd.get()));
            }
        }

        // Write all L1X metadata.
        for (size_t il1x=0 ; il1x<il1x_start ; il1x++) {
            handle(l1x[il1x]->write_metadata(
                frm.get() // Frame to write L1X from.
            ));
        }
        handle(l1a->write_metadata(
            frm.get() // Frame to write L1X from.
        ));

        if (il1x_end > L1X_RAD) continue;
        if (il1x_start > L1X_RAD) {
            // Step: Modulate. interpolated truth to L1X_RAD.
            // L1X output RAD.
            handle(l1x[L1X_RAD]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));
        }

        if (il1x_end > L1X_FOV) continue;
        if (il1x_start > L1X_FOV) {

            // Step: Uncalibate spectra. L1X_RAD to L1X_FOV.
            handle(frm->uncalibrate_spectra(
                ckd.get()
            ));
            // L1X output FOV.
            handle(l1x[L1X_FOV]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        if (il1x_end > L1X_REBIN) continue;
        if (il1x_start > L1X_REBIN) {

            // Step: Draw the spectra on the detector.
            handle(frm->draw_on_detector(ckd.get()));

            // L1X output REBIN.
            handle(l1x[L1X_REBIN]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

         }

        if (il1x_end > L1X_STRAY) continue;
        if (il1x_start > L1X_STRAY) {
            //added there :
            // Step: Straylight.
            //handle(frm->apply_straylight(ckd.get()));
            //end of the additions
            // L1X output STRAY.
            handle(l1x[L1X_STRAY]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        if (il1x_end > L1X_UNBIN) continue;
        if (il1x_start > L1X_UNBIN) {

            // Step: Straylight.
            handle(frm->apply_straylight(set_main, ckd.get()));

            // L1X output UNBIN.
            handle(l1x[L1X_UNBIN]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));
        }

        if (il1x_end > L1X_PRNU) continue;
        if (il1x_start > L1X_PRNU) {

            // L1X output PRNU.
            handle(l1x[L1X_PRNU]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        if (il1x_end > L1X_NONLIN) continue;
        if (il1x_start > L1X_NONLIN) {

            // step: prnu.
            handle(frm->apply_prnu(ckd.get()));

            // L1X output NONLIN.
            handle(l1x[L1X_NONLIN]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        if (il1x_end > L1X_NOISE) continue;
        if (il1x_start > L1X_NOISE) {

            // Step: Non-linearity.
            handle(frm->apply_nonlinearity(ckd.get()));

            // L1X output NOISE.
            handle(l1x[L1X_NOISE]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        if (il1x_end > L1X_DARK) continue;
        if (il1x_start > L1X_DARK) {

            // Step: Dark current.
            handle(frm->apply_dark_current(ckd.get()));

            // L1X output DARK.
            handle(l1x[L1X_DARK]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));

        }

        // Compute and print the optimal coadding factor and exposure time
        if (set_main.t_dwell > 0.0) {
            double max_val {};
            for (int i {}; i < ckd->npix; ++i) {
                if (ckd->dark_current[i] < -100.0
                    || ckd->dark_current[i] > 550.0) {
                // if (!ckd->mask[i]) {
                    max_val = std::max(max_val, frm->image[i]);
                }
            }
            const double I_sig { max_val / frm->exposure_time };
            const double coadd_raw {
                I_sig * set_main.t_dwell / (set_main.f_sat * set_main.full_well
                                            + I_sig * set_main.t_dead)
            };
            const int n_coadd { static_cast<int>(std::ceil(coadd_raw)) };
            std::cout << "Optimal coadding factor (float) " << coadd_raw
                      << '\n';
            std::cout << "Optimal coadding factor (int) " << n_coadd << '\n';
            std::cout << "Optimal exposure time "
                      << (set_main.t_dwell / n_coadd) - set_main.t_dead
                      << " s\n";
        }
        if (il1x_end > L1X_RAW) continue;
        if (il1x_start > L1X_RAW) {

            // Step: Dark offset.
            handle(frm->apply_dark_offset(ckd.get()));

            // L1X output RAW.
            handle(l1x[L1X_RAW]->write(
                frm.get(), // Frame to write L1X from.
                &set_ns, // Noise settings.
                ckd.get() // Calibration key data (for detector and spectrum dimensions).
            ));
        }

        if (il1x_end > L1X_L1A) continue;
        // if-clase il1x_start greater then L1X_L1A. That is always true,
        // because you cannot start with L1A.

        // Step: Add noise.
        if (set_ns.noise_apply) {
            static std::mt19937 gen { set_ns.seed };
            for (int i {}; i < ckd->npix; ++i) {
                const double noise_value {
                    std::sqrt(ckd->noise_n[i] * ckd->noise_n[i]
                              + frm->image_with_current[i] * ckd->noise_g[i])
                };
                std::normal_distribution<> d { 0.0, noise_value };
                frm->image[i] += d(gen);
            }
        }

        // Step: Perform the binning.
        if (set_main.binning_table_id != 0) {
            vector<double> image_binned(bin->npix,0);
            for (size_t ipix_unbinned=0 ; ipix_unbinned<ckd->npix ; ipix_unbinned++) {
                uint32_t &ipix = bin->pixelpointer[ipix_unbinned];
                if (ipix < bin->npix) {
                    image_binned[ipix] += frm->image[ipix_unbinned];
                }
            }
            frm->image = image_binned;
        }

        // Step: Multiply with co-additions.
        for (size_t ipix=0 ; ipix<npix_binned ; ipix++) {
            frm->image[ipix] *= frm->nr_coadditions;
        }

        // Step: Convert to unsigned shorts.
        frm->image_ints.resize(npix_binned);
        for (size_t ipix=0 ; ipix<npix_binned ; ipix++) {
            frm->image_ints[ipix] = static_cast<uint16_t>(frm->image[ipix]+0.5);
        }

        // L1X output L1A.
        handle(l1a->write(
            frm.get(), // Frame to write L1A from.
            &set_ns, // Noise settings (not used for L1A).
            ckd.get() // Calibration key data (not used for L1A).
        ));

    }

    return 0;

} // }}}

