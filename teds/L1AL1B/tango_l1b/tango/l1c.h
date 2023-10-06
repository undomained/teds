// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef L1C_H
#define L1C_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_l1c : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_l1c(
        Logger *creator
    );
    // Specific settings.
    string demfile = ""; // Detailed elevation map. Leave empty for placeholder DEM function.
    double semi_major_axis = NC_FILL_DOUBLE; // Earth property (move to database).
    double semi_minor_axis = NC_FILL_DOUBLE; // Earth property (move to database).
    double latitude_tol = NC_FILL_DOUBLE; // Tolerance for iterative calculation of latitude during geolocation (radians).
    string l1bfile = ""; // Input L1B file (no L1A files are used here).
    string l1cfile = ""; // Output L1C file.
    vector<double> lat; // Latitudes of target positions.
    vector<double> lon; // Longitudes of target positions.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class L1C : public Processor { // {{{

    // A public constructor.
    public:
    L1C(
        Logger *creator,
        CKD *ckd_arg
    );

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_l1c> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;

}; // }}}

#endif
