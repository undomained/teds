#include "header.h"
#include "functions.h"
#include "logger.h"
#include "program.h"

// Argument switches. {{{
const size_t nswitch = 3;
const size_t iswitch_verbose = 0; // -v
const size_t iswitch_foldsettings = 1; // -f
const size_t iswitch_timestamps = 2; // -t
const vector<string> switches = {
    "-v", // Verbose.
    "-f", // Fold settings.
    "-t" // Time staps.
};
const vector<string> switch_meanings = {
    "More logging on screen.", // Verbose (-v).
    "Add vim fold markers around settings in log file.", // Fold settings (-f).
    "Add time stamps in file names of any output file." // Time staps (-t).
};
// }}}

void printusage( // {{{
    const char *executable
)
{
    string programname = executable;
    size_t pos = programname.find("/");
    while (pos != string::npos) {
        programname = programname.substr(pos+1);
        pos = programname.find("/");
    }
    string all_switches = "";
    for (size_t iswitch=0 ; iswitch<nswitch ; iswitch++) {
        all_switches += format(" [%s]",switches[iswitch].c_str());
    }
    printf("%s%s settings_file.cfg\n",programname.c_str(),all_switches.c_str());
    for (size_t iswitch=0 ; iswitch<nswitch ; iswitch++) {
        printf("  [%s]   %s\n",switches[iswitch].c_str(),switch_meanings[iswitch].c_str());
    }
} // }}}

int main( // {{{
    int argc,
    char *argv[]
)
{
    vector<bool> flags(nswitch,false);

    // Regular argument.
    string settings_file = "";

    for (int iarg=1 ; iarg<argc ; iarg++) {
        bool switch_found = false;
        for (size_t iswitch=0 ; iswitch<nswitch ; iswitch++) {
            if (switches[iswitch].compare(argv[iarg]) == 0) {
                if (flags[iswitch]) {
                    printusage(argv[0]);
                    return 1;
                }
                flags[iswitch] = true;
                switch_found = true;
            }
        }
        if (!switch_found) {
            if (settings_file.compare("") != 0) {
                printusage(argv[0]);
                return 1;
            }
            settings_file = argv[iarg];
        }
    }
    if (settings_file.compare("") == 0) {
        printusage(argv[0]);
        return 1;
    }

    // Create time stamp.
    string timestamp = "";
    if (flags[iswitch_timestamps]) {
        // This is a slightly adapted version of the function now_timestring.
        // It has a different format and therefore also a different length.
        time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        struct tm *timeinfo = localtime(&now);
        const size_t sz = 17; // With \0.
        char buffer[sz];
        strftime(buffer,sz,"_%Y%m%dT%H%M%S",timeinfo);
        timestamp = string(buffer);
    }

    // Derive log file.
    size_t dotpos = settings_file.find(".");
    string log_file = settings_file.substr(0,dotpos) + timestamp + ".txt";
    if (settings_file.compare(log_file) == 0) {
        printf("Error: Settings file should never have '.txt' extension, because that interferes with the log file, unless you switch on time stamps.\n");
        return 1;
    }

    Logger trunk(log_file,flags[iswitch_verbose],timestamp);
    Program prog(&trunk);

    // Execute the program.
    return prog.execute(settings_file,flags[iswitch_foldsettings]);

} // }}}

