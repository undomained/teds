#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

class UTC : public Logger {

    // This should have been done with a search function and private members.
    public:
    UTC(
        Logger *creator
    );
    ~UTC();

    size_t dim_utc;
    vector<int> mjd_utc;
    vector<double> tdiff_utc;

    int init(
        string &utcfile
    );
};

} // namespace tango
