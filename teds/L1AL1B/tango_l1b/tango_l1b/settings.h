// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Defines an abstract class for storing all user defined
// configuration parameters. It is initialized using a YAML
// configuration file. The parameter values are defined in a derived
// class. Each configuration parameter is a Setting instance for which
// we define its YAML path, default value (omitted for
// Setting<std::optional<T>>), an info string, and optionally whether
// it's a normal or debug variable (level).
//
// class SettingsDerived : public Settings
// {
// public:
//     SettingsDerived(const std::string& yaml_file) : Settings { yaml_file } {}
//     struct
//     {
//         Setting<bool> run { { "main", "run" }, true, "description" };
//     } main;
//     auto scanKeys() -> void override
//     {
//         scan(main.run);
//     }
//     auto checkParameters() -> void override {}
// };
//
// Storing the parameters in (nested) structs is a matter of
// organization and not required.
//
// The level shows whether the parameter comes with a safe default
// value that should only be changed for debugging/experimental
// purposes (debug) or whether it's a normal parameter requiring user
// attention (normal).
//
// Notes
//
// Any newline symbols in the Setting info string are recognized by
// YAML::Emitter and are good for readability when printing the
// configuration.

#pragma once

#include "yaml.h"

namespace tango {

class Settings
{
private:
    // Whether to load or dump (for output) a YAML node in the scan function
    bool do_dump { false };
    // Check for unrecognized keywords
    auto unrecognizedKeywordCheck() const -> void;
    // Keeps track of the location of the current map which equals the
    // YAML keys of a setting minus the last key. Example: for the
    // setting pol.print_debug.enabled the map location is { "pol",
    // "print_debug" }. If there is a change in cur_map_loc then a
    // YAML::BeginMap or YAML::EndMap must be inserted to
    // YAML::Emitter.
    std::vector<std::string> cur_map_loc {};
    // An Emitter instance for converting the configuration into a string
    Emitter yaml_emitter {};
    // Feed a C++ data structure (a Setting instance) into
    // YAML::Emitter which converts it into a character
    // stream. automatically insert YAML::BeginMap and YAML::EndMap
    // where necessary.
    template <typename T>
    auto dump(Emitter& emitter, const Setting<T>& setting) -> void
    {
        const auto& keys { setting.yaml_keys };
        // We define the YAML keys minus the last key of a setting as
        // its location. If this is shorter than cur_map_loc then we
        // are at the end of map(s) and one or more YAML::EndMap must
        // be inserted.
        while (cur_map_loc.size() + 1 > keys.size()) {
            emitter << YAML::EndMap;
            cur_map_loc.pop_back();
        }
        // If the current map location is shorter or equal to the
        // setting location, compare the keys starting with the last
        // item. If they differ, the map has changed and one or more
        // YAML::EndMap must be inserted.
        for (int i { static_cast<int>(
               std::min(cur_map_loc.size(), keys.size() - 1) - 1) };
             i >= 0;
             --i) {
            if (cur_map_loc.at(i) != keys.at(i)) {
                emitter << YAML::EndMap;
                cur_map_loc.pop_back();
            }
        }
        // At this point, up to the size of cur_map_loc, the current
        // map location agrees with the setting location. If the
        // setting location is longer then we are starting a new map
        // and the YAML::Key and YAML::BeginMap directives must be
        // inserted accordingly.
        for (int i { static_cast<int>(cur_map_loc.size()) };
             i < static_cast<int>(keys.size() - 1);
             ++i) {
            cur_map_loc.push_back(keys.at(i));
            emitter << YAML::Key << cur_map_loc.back() << YAML::Value
                    << YAML::BeginMap;
        }
        // Finally, we have sorted out the correct nesting of maps
        // (you can think of this as the correct indentation in a YAML
        // file) and can print the value (which can consist of several
        // details) of the setting.
        emitter << setting;
    }

protected:
    // Full configuration as read from file
    YAML::Node config {};
    // Full configuration with default values
    YAML::Node default_config {};
    // The scan function records every key it comes across in
    // all_valid_keys. This can be useful later for recognizing
    // unknown keys.
    std::vector<std::vector<std::string>> all_valid_keys {};
    // Search for invalid values of configuration parameters
    // (e.g. whether a location points to a real file or whether some
    // number is within a reasonable range) and inconsistencies
    // between parameters. Editing parameter values is allowed and
    // thus will not be declared const here.
    virtual auto checkParameters() -> void = 0;

public:
    Settings() = default;
    Settings(const std::string& yaml_file)
      : config { YAML::LoadFile(yaml_file) }
    {}
    Settings(const Settings& /* settings */) {};
    // Read input from a YAML configuration file and check parameters
    // for correctness.
    auto init() -> void;
    // The main function for processing YAML input. This is called in init.
    virtual auto scanKeys() -> void = 0;
    // Find a YAML node and either dump into a string or convert it
    // into a C++ data structure
    template <typename T>
    auto scan(Setting<T>& item)
    {
        if (do_dump) {
            dump(yaml_emitter, item);
            return;
        }
        if (item.yaml_keys.empty()) {
            return;
        }
        all_valid_keys.push_back(item.yaml_keys);
        YAML::Node node { YAML::Clone(config) };
        // Search for a node in the YAML configuration. If found, set
        // item to that value. Otherwise, return and leave the default
        // value unmodified.
        for (const auto& key : item.yaml_keys) {
            node = node[key];
            if (!node) {
                return;
            }
        }
        try {
            item = node.as<T>();
        } catch (const YAML::BadConversion&) {
            // Let the user know exactly which key failed to convert
            std::string str_value {};
            try {
                str_value = node.as<std::string>();
            } catch (const YAML::BadConversion&) {
                // If can't extract the value then the exception
                // message will be less informative.
            }
            throw std::runtime_error { "cannot set " + item.keyToStr()
                                       + ", which is of type " + item.type
                                       + ", to the value " + str_value };
        }
    }
    // Convert current configuration into a string. If only the
    // default constructor was called, print the default
    // configuration.
    auto c_str(const bool verbose = true) -> const char*;
    // Return a clone of the configuration. If a configuration
    // variable was not set by the user, the default value is shown
    // instead.
    auto getConfig() const -> std::string;
    virtual ~Settings() = default;
};

} // namespace tango
