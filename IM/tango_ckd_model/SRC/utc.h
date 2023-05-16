#ifndef UTC_H
#define UTC_H

#include "header.h"
#include "logger.h"

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

#endif
