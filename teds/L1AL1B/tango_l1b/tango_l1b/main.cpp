// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "driver.h"

// Argument switches. {{{
const size_t nswitch = 3;
const size_t iswitch_verbose = 0; // -v
const size_t iswitch_foldsettings = 1; // -f
const size_t iswitch_timestamps = 2; // -t
const std::vector<std::string> switches = {
    "-v", // Verbose.
    "-f", // Fold settings.
    "-t" // Time stamps.
};
const std::vector<std::string> switch_meanings = {
    "More logging on screen.", // Verbose (-v).
    "Add vim fold markers around settings in log file.", // Fold settings (-f).
    "Add time stamps in file names of any output file." // Time stamps (-t).
};
// }}}

void printusage( // {{{
    const char *executable
)
{
    std::string programname = executable;
    size_t pos = programname.find("/");
    while (pos != std::string::npos) {
        programname = programname.substr(pos+1);
        pos = programname.find("/");
    }
    std::string all_switches = "";
    for (size_t iswitch=0 ; iswitch<nswitch ; iswitch++) {
        all_switches += std::format(" [%s]",switches[iswitch].c_str());
    }
    printf("%s%s settings_file.cfg\n",programname.c_str(),all_switches.c_str());
    for (size_t iswitch=0 ; iswitch<nswitch ; iswitch++) {
        printf("  [%s]   %s\n",switches[iswitch].c_str(),switch_meanings[iswitch].c_str());
    }
} // }}}

std::string read_key_from_config(// {{{
    std::string settings_file,
    std::string key
)
{
    std::ifstream fileStream(settings_file.c_str());
    std::string line;
    while (getline(fileStream, line))
        {
            std::istringstream is_line(line);
            std::string found_key;
            if (std::getline(is_line, found_key, '='))
            {
                // Remove trailing spaces from key.
                while (found_key[found_key.size()-1] == ' ') found_key.erase(found_key.size()-1,1);
                // Remove leading spaces from key.
                while (found_key[0] == ' ') found_key.erase(0,1);
                std::string value;
                if (found_key[0] == '#')
                    // Comment line. Skip it
                    continue;
                if (strcmp(found_key.c_str(), key.c_str()) == 0){
                    // Found the key we are looking for
                    if (std::getline(is_line, value))
                    {
                        // Remove leading spaces from value.
                        while (value[0] == ' ') value.erase(0,1);
                        // Remove trailing spaces from value.
                        while (value[value.size()-1] == ' ') value.erase(value.size()-1,1);
                        return value;
                    }
                }
            }
        }
    printf("Warning key '%s' not found in settings file (%s)! Returning empy string!\n", key.c_str(), settings_file.c_str());
    return "";
} // }}}

std::string get_log_file_name(// {{{
    std::string settings_file, 
    std::string log_file_path,
    std::string timestamp
)
{

    // Derive log file name.
    // Name of logfile is derived from name the settings file but with extension .log instead of .cfg
    // log file path is taken from configuration. If it is not set there the log file will be written to
    // the directory from which job is run (should be teds directory)

    // Get position of the last / in the name of the settings file
    size_t slashpos = settings_file.rfind("/");
    // Get position of the start of the file name (i.e. without the path)
    size_t file_name_start = 0;
    if (slashpos != std::string::npos) {
        file_name_start = slashpos+1;
    }
    // Get the file name
    std::string file_name = settings_file.substr(file_name_start);

    // Look for the extension dot.
    size_t dotpos = file_name.find(".");
    // Give warning and exit if settings file is without extension (then you are not running the code as described in Readme).
    if (dotpos == std::string::npos) {
        printf("Warning: Settings file has no extension. Please provide a settings file WITH extension .cfg!\n");
        return "";
    }

    std::string log_file = log_file_path + file_name.substr(0,dotpos) + timestamp + ".log";
    printf("logfile name: {%s}\n", log_file.c_str());
    return log_file;

} // }}}



int main( // {{{
    int argc,
    char *argv[]
)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    std::vector<bool> flags(nswitch,false);

    // Regular argument.
    std::string settings_file = "";

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
    std::string timestamp = "";
    if (flags[iswitch_timestamps]) {
        // This is a slightly adapted version of the function now_timestring.
        // It has a different format and therefore also a different length.
        time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        struct tm *timeinfo = localtime(&now);
        const size_t sz = 17; // With \0.
        char buffer[sz];
        strftime(buffer,sz,"_%Y%m%dT%H%M%S",timeinfo);
        timestamp = std::string(buffer);
    }

    // At this moment no logger defined so not possible to use settings objects
    // Use fct to find log_file_path in config file
    std::string log_file_path = read_key_from_config(settings_file, "log_file_path");

    std::string log_file = get_log_file_name(settings_file, log_file_path, timestamp);
    if (log_file.empty()) return 1;

//    // Derive log file.
//    // The logfile is equal to the settings file, except that the extension
//    // is substituted to 'txt'. In principle, the settings file should be
//    // <filename>.cfg, but it may be possible that it has a different
//    // extension or no extension at all. In the latter case, we just
//    // add `.txt' to the extensionless filename for the logfile.
//    size_t dotpos = settings_file.rfind(".");
//    // It is possible that the settings file has no extension and that
//    // a relative path with .. is included. Than, there is a rightmost dot
//    // that should not be interpreted as the start of an extension. To
//    // prevent that, we ignore the found position if it is left of the
//    // rightmost slash.
//    size_t slashpos = settings_file.rfind("/");
//    if (slashpos != string::npos && dotpos < slashpos) dotpos = string::npos;
//    string log_file = settings_file.substr(0,dotpos) + timestamp + ".txt";
//    if (settings_file.compare(log_file) == 0) {
//        printf("Error: Settings file should never have '.txt' extension, because that interferes with the log file, unless you switch on time stamps.\n");
//        return 1;
//    }

    tango::Logger trunk(log_file,flags[iswitch_verbose],timestamp);
    int stat;
    std::string instrument_cal_choice = read_key_from_config(settings_file, "instrument_cal");
    printf("Info: instrument choice: --%s-- \n",instrument_cal_choice.c_str());
    if (strcmp(instrument_cal_choice.c_str(), "spexone") == 0){
        tango::Spexone_cal prog(&trunk);
        // Execute the program.
        stat = prog.execute(settings_file,flags[iswitch_foldsettings]);
    }


#ifdef USE_MPI
    if (stat) {
        MPI_Abort(MPI_COMM_WORLD, stat);
    } else {
        MPI_Finalize();
    }
#endif
    return stat;
} // }}}
