// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "batch.h"

Batch::Batch() {}
Batch::~Batch() {}

// Batch resizer.
void Batch::setSize( // {{{
    size_t a_nl1a
)
{
    nl1a = a_nl1a;
    l1a.resize(nl1a);
    remove.resize(nl1a,true); // Default remove flag.
} // }}}

