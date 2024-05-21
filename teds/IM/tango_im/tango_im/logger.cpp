#include "header.h"
#include "functions.h"
#include "logger.h"

namespace tango {

Logger::Logger( // {{{
    Logger *creator
)
{
    // Take over the necessary information from the creator.
    // This preserves eventual changes, most likely in the name.
    verboseflag = creator->verboseflag; // Can be adapted.
    filename = creator->filename; // I do not think you want to use this.
    file = creator->file; // This copies a pointer.
    name = creator->name; // Can be adapted.
    nwarn = 0; // Own warnings first.
    dead = false; // Most loggers are born alive.
    perc = creator->perc; // Copy shared pointer.
    timestamp = creator->timestamp; // Copy time stamp for output NetCDF.
    // But we need the pointer to the trunk.
    if (creator->trunk == NULL) {
        // Your creator is the real trunk
        trunk = creator;
    } else {
        // Take the same trunk as the creator.
        trunk = creator->trunk;
    }
} // }}}
Logger::Logger( // {{{
    string a_filename,
    bool a_verboseflag,
    string a_timestamp
)
{

    trunk = NULL; // We are the trunk.
    verboseflag = a_verboseflag;
    filename = a_filename;
    file = fopen(filename.c_str(),"w");
    name = "main";
    nwarn = 0;
    dead = false;
    perc = make_shared<vector<Percentage_struct>>();
    timestamp = a_timestamp;
    writelog(log_verbose,"Opened log file: '%s'.",filename.c_str());

} // }}}
Logger::~Logger( // {{{
)
{
    if (trunk == NULL) {

        // Deconstruct the logger.
        if (dead) {
            if (nwarn == 1) writelog(log_info,"There was 1 warning.");
            else if (nwarn > 1)  writelog(log_info,"There were %zu warnings.",nwarn);
            writelog(log_info,"Program exited with error.");
        } else {
            if (nwarn == 0) writelog(log_info,"Program successfully finished.");
            else if (nwarn == 1) writelog(log_info,"Program finished with 1 warning.");
            else writelog(log_info,"Program finished with %zu warnings.",nwarn);
        }
        writelog(log_verbose,"Closing log file: '%s'.",filename.c_str());
        fclose(file);

    } else {

        // Transfer information to trunk.
        trunk->nwarn += nwarn;
        if (dead) trunk->dead = true; // Or-construction.

        // Do not close off the file.

    }

} // }}}
void Logger::setName( // {{{
    string a_name
)
{
    name = a_name;
} // }}}
void Logger::writelog( // {{{
    loglevel_t level,
    string fmt,
    ...
)
{

    if (level == log_verbose && !verboseflag) return;

    va_list ap;
    va_start(ap,fmt);
    va_list ap_cpy;
    va_copy(ap_cpy,ap);
    size_t sz = vsnprintf(NULL,0,fmt.c_str(),ap) + 1;
    string msg(sz,' ');
    vsnprintf(&msg.front(),sz,fmt.c_str(),ap_cpy);
    va_end(ap_cpy);
    // Test with name.
    string cat;
    if (level == log_verbose) cat = "V";
    if (level == log_trace) cat = "T";
    if (level == log_debug) cat = "D";
    if (level == log_info) cat = "I";
    if (level == log_warning) cat = "W";
    if (level == log_error) cat = "E";
    if (level == log_verbose || level == log_trace || level == log_info || level == log_warning || level == log_error) {
        size_t nperc = perc->size();
        if (nperc != 0) printf("\n");
        printf("[%s] %s: %s\n",cat.c_str(),name.c_str(),msg.c_str());
        if (nperc != 0 && level != log_error) {
            Percentage_struct &perc_cur = (*perc)[nperc-1];
            printf("%s: %d%%",perc_cur.percentage_message.c_str(),perc_cur.percentage);
            cout.flush();
        }
    }
    if (level == log_debug || level == log_info || level == log_warning || level == log_error) {
        fprintf(file,"[%s] %s: %s\n",cat.c_str(),name.c_str(),msg.c_str());
        fflush(file);
    }
    if (level == log_warning) nwarn++;
    if (level == log_error) {
        dead = true;
        perc->resize(0); // Remove any percentages when logging an error.
    }
    va_end(ap);
} // }}}

void Logger::percentagelog_open( // {{{
    string fmt,
    ...
)
{

    // The differences with a normal log.
    // The log level is always log_trace.
    // And the \n is omitted and thus a flush is called.
    // A ": 0%" is added.
    // The percentage flag is set.
    // If the percentage is already set, it is an error.
    va_list ap;
    va_start(ap,fmt);
    va_list ap_cpy;
    va_copy(ap_cpy,ap);
    size_t sz = vsnprintf(NULL,0,fmt.c_str(),ap) + 1;
    string msg(sz,' ');
    vsnprintf(&msg.front(),sz,fmt.c_str(),ap_cpy);
    va_end(ap_cpy);
    size_t nperc = perc->size();
    if (nperc != 0) printf("\n");
    // Append the new percentage structure.
    perc->resize(nperc+1);
    Percentage_struct &perc_cur = (*perc)[nperc];
    // Saves message and zero progress.
    perc_cur.percentage = 0;
    perc_cur.percentage_message = format("[T] %s: %s",name.c_str(),msg.c_str());
    printf("%s: %d%%",perc_cur.percentage_message.c_str(),perc_cur.percentage);
    cout.flush();
    va_end(ap);
} // }}}
void Logger::percentagelog_progress( // {{{
    size_t idx,
    size_t sz
)
{
    size_t nperc = perc->size();
    if (nperc == 0) return;
    Percentage_struct &perc_cur = (*perc)[nperc-1];
    // Just update the percentage.
    int percentage_new = 100*idx/sz;
    if (percentage_new > perc_cur.percentage) {
        if (perc_cur.percentage > 9) printf("\b");
        printf("\b\b%d%%",percentage_new); cout.flush();
        perc_cur.percentage = percentage_new;
    }
} // }}}
void Logger::percentagelog_close( // {{{
    int percentage_end
)
{
    size_t nperc = perc->size();
    Percentage_struct &perc_cur = (*perc)[nperc-1];
    if (perc_cur.percentage > 9) printf("\b");
    printf("\b\b%d%%\n",percentage_end);
    perc->resize(nperc-1);
    if (nperc == 1) return;
    Percentage_struct &perc_prv = (*perc)[nperc-2];
    printf("%s: %d%%",perc_prv.percentage_message.c_str(),perc_prv.percentage);
    cout.flush();
} // }}}

} // namespace tango
