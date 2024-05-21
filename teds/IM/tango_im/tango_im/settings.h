#ifndef SETTINGS_H
#define SETTINGS_H

#include "header.h"
#include "logger.h"

// Settings recognition macros. {{{

// These also use macros because it would be called with a lot of same arguments.
// And it is no rocket science. It is just an instruction to recognize a settings.
// Moreover, it uses the hash tag to convert a variable name to a string and that
// is the system that we want to use.

// Recognize regular setting.
#define recognize_setting(option) \
    if (key == #option) \
    { \
        stream << value; \
        stream >> option; \
        recognized = true; \
    }
// A vector.
#define recognize_setting_vector(option) \
    { \
        if (key.find(".") == string::npos) { \
            if (key == #option) { \
                stream << value; \
                size_t el = option.size(); \
                option.resize(el+1); \
                stream >> option[el]; \
                recognized = true; \
            } \
        } else { \
            istringstream keystream(key); \
            string vecstring; \
            if (getline(keystream,vecstring,'.')) { \
                if (vecstring == #option) { \
                    string elstring; \
                    stringstream elstream; \
                    if (getline(keystream,elstring)) { \
                        size_t el; \
                        elstream << elstring; \
                        elstream >> el; \
                        stream << value; \
                        if (option.size() < el+1) option.resize(el+1); \
                        stream >> option[el]; \
                        recognized = true; \
                    } \
                } \
            } \
        } \
    }

// }}}

// Settings class. These are one class for processor settings and main settings.
class Settings : public Logger { // {{{

    // Protected constructor to prevent any raw Settings instances to be created.
    protected:
    Settings(
        Logger *creator
    );

    public:
    // Virtual destructor to prevent surprises.
    // If a subclass object is destructed when being a Settings pointer, the
    // destructor of the subclass is only called when the Settings object has a
    // virtual destructor. Note that the virtual destructor itself is also called
    // when a subclass has its own destructor.
    virtual ~Settings() {}

    // This routine opens the settings file and reads contents. A virtual routine
    // is called to add the calibration-step specific settings.
    int init(
        string &settings_file // Settings file name.
    );

    string tag; // Tag that must be searched when reading the settings. Will be set by legal child or grandchild classes.

    // Routine to overwrite, where either common processor or logistic (main) settings are read.
    protected:
    virtual int init_common(
        stringstream &stream,
        string &key,
        string &value,
        bool &recognized
    ) = 0;

}; // }}}

#endif
