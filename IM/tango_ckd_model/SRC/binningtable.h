#ifndef BINNINGTABLE_H
#define BINNINGTABLE_H

#include "header.h"
#include "logger.h"

// Forward declaration.
class NetCDF_object;

// Binning table file.
class Binningtable_file : public Logger { // {{{

    public:
    size_t dim_binningtable_spat;
    size_t dim_binningtable_spec;
    size_t prebinning_spat;
    size_t prebinning_spec;
    size_t dim_detector_spat;
    size_t dim_detector_spec;
    unique_ptr<NetCDF_object> nc;

    Binningtable_file(
        Logger *creator
    );
    ~Binningtable_file();

    int read(
        string &binningtable_filename, // File with all the used binning tables.
        size_t a_dim_detector_spat,
        size_t a_dim_detector_spec
    );

}; // }}}

// Binning table.
class Binningtable : public Logger { // {{{

    public:
    Binningtable(
        Logger *creator
    );
    ~Binningtable();

    vector<uint32_t> pixelpointer; // The actual binning table. Not allocated for no binning (if trivial is true).

    // Read the binning table from a NetCDF file.
    int read(
        Binningtable_file *bin_file, // Binning table file.
        uint8_t id // Mask identifier.
    );

    size_t npix_unbinned; // Because used often in own routines.
    // Number of pixels: This can be less than the L1A number of pixels in which case
    // the last pixels need to be skipped.
    size_t npix; // Also used often.

}; // }}}

#endif
