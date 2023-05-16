// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef BINNINGTABLE_H
#define BINNINGTABLE_H

#include "header.h"
#include "logger.h"

// Forward declaration.
class NetCDF_object;
class CKD;

// Binning table file.
class Binningtable_file : public Logger { // {{{

    public:
    size_t dim_binningtable_spat;
    size_t dim_binningtable_spec;
    size_t prebinning_spat;
    size_t prebinning_spec;
    unique_ptr<NetCDF_object> nc;

    Binningtable_file(
        Logger *creator
    );
    ~Binningtable_file();

    int read(
        string &binningtable_filename, // File with all the used binning tables.
        CKD *ckd // Pointer to the general CKD.
    );

}; // }}}

// Binning table.
class Binningtable : public Logger { // {{{

    public:
    Binningtable(
        Logger *creator
    );
    ~Binningtable();

    bool trivial; // Flag for trivial (no binning) table.
    vector<uint32_t> pixelpointer; // The actual binning table. Not allocated for no binning (if trivial is true).
    vector<uint32_t> binsizes; // Number of small pixels included in each bin. Not allocated for no binning (if trivial is true).

    // The used binned CKD, will point to own_binned_ckd for non-trivial
    // binning tables and to the general CKD for trivial binning tables.
    CKD *binned_ckd;

    // Getter.
    size_t getBinsize(
        size_t ipix
    );

    // Read the binning table from a NetCDF file.
    int read(
        Binningtable_file *bin_file, // Binning table file.
        uint8_t id, // Mask identifier.
        CKD *ckd // Pointer to the general CKD.
    );

    size_t npix_unbinned; // Because used often in own routines.
    // Number of pixels: This can be less than the L1A number of pixels in which case
    // the last pixels need to be skipped.
    size_t npix; // Also used often.

    private:
    int average_ckd(
        size_t dim_slower, // Product of dimensions slower than pixels.
        size_t dim_quicker, // Product of dimensions quicker than pixels.
        vector<double> &src, // Source vector from original CKD.
        vector<double> &dest // Destination vector from binned CKD.
    );

    int explore_cell(
        size_t ispat,
        size_t ispec_start,
        size_t &ispec_end,
        size_t &ipix,
        size_t &ispat_start,
        size_t &ispat_end
    );

    // Binned calibration key data (only for detector steps).
    unique_ptr<CKD> own_binned_ckd; // Not used for trivial.
    // Public pointer will point to this one or the original CKD.

}; // }}}

#endif
