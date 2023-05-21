#include "header.h"
#include "functions.h"
#include "logger.h"
#include "settings.h"

// Constructor, only constructing the logger part.
Settings::Settings( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}

// Read the settings file and stores the contents in the structure.
int Settings::init( // {{{
    string &settings_file // File name for the settings.
)
{

    // Turn file into a stream.
    ifstream setstream(settings_file.c_str());
    check_error(setstream.fail(),"Error opening settings file: '%s'.",settings_file.c_str());

    // First, we will browse to the correct tag. Then, we will continue until the
    // end of the file or until the next tag comes.
    bool tag_found = false;
    string tag_brackets = "[" + tag + "]";

    // Loop per line.
    string line;
    while (getline(setstream,line)) {

        // Clip off character 13 at the end of a line.
        if (line.length() > 0 && line[line.length()-1] == 13) line = line.substr(0,line.length()-1);

        // Clip off any comments. Comments are prefixed with //, # or ;.
        // So, it is Python, C++ or IDL style. At least, that is what the
        // default Vim syntax highlighting rules for cfg-files suggest.
        size_t found = line.find("#");
        line = line.substr(0,found);
        found = line.find("//");
        line = line.substr(0,found);
        found = line.find(";");
        line = line.substr(0,found);
        // Remove leading spaces.
        while (line.length() > 0 && line[0] == ' ') line = line.substr(1);
        // Remove trailing spaces.
        while (line.length() > 0 && line[line.length()-1] == ' ') line = line.substr(0,line.length()-1);
        // Ignore if it is empty.
        if (line.length() == 0) continue;

        // Browse to tag.
        if (!tag_found) {
            if (line.compare(tag_brackets) == 0) tag_found = true;
            continue;
        }
        // Stop when new tag is detected.
        if (line[0] == '[') break; // New tag found, ending ] not checked.

        istringstream linestream(line);
        string key;
        // This if-clause is true if there is an equal sign. The key will be filled
        // with everything left of the equal sign.
        if (getline(linestream,key,'=')) {

            // Remove trailing spaces from key.
            while (key[key.size()-1] == ' ') key.erase(key.size()-1,1);

            string value; // String after the equal sign, to be interpreted to the right type.
            stringstream stream; // Temporarily-used string.

            // Read the rest of the line and write down the value. Since we already read
            // up to the equal sign, this is the part right of the equal sign.

            if (getline(linestream,value)) {

                // Remove leading spaces from value.
                while (value[0] == ' ') value.erase(0,1);

                // Flag for having recognized the setting. If nothing is recognized, it
                // is a wrong line, possibly a misspelled key.
                bool recognized = false;

                // Read processor settings or main settings.
                handle(init_common(stream,key,value,recognized));

                // Handle error for no recognition.
                check_error(!recognized,"Settings error: Invalid key: %s",key.c_str());

            } else {
                // No value found.
                raise_error("Settings error. No proper value given after key: %s",key.c_str());
            }
        } else {
            // No key found. That is, no equal sign.
            raise_error("Settings error. No '=' found in line: %s",line.c_str());
        }
    }

    check_error(!tag_found,"Error: tag [%s] not found in settings file.",tag.c_str());

    return 0;

} // }}}
