// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef BATCH_H
#define BATCH_H

#include "header.h"

// Forward declaration.
class L1A;

// Organization of a batch of L1A images that are processed at the same time.
struct Batch { // {{{

    // Do-nothing constructor and destructor.
    Batch();
    ~Batch();

    // All these vectors have the same number of elements.
    size_t nl1a; // Size of all these vectors.
    vector<L1A *> l1a; // Pointers to L1A instances.
    size_t ipix_start = 0; // Where the fraction starts.
    size_t ipix_end = NC_FILL_UINT64; // Where the fraction ends (fillvalue will be interpreted as l1a->bin->npix).
    // Recycling is possible if the same image is needed more than once.
    // Only if the second fraction is a subset of the first fraction, this is possible.
    // When no co-adding is requested, recycling is only possible if the fractions are
    // identical.
    vector<bool> remove; // Flag to remove image after use.
    void setSize(
        size_t nl1a
    );

}; // }}}

#endif
