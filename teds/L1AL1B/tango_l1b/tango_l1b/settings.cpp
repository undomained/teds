// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "settings.h"

#include <iostream>

namespace tango {

auto Settings::init() -> void
{
    default_config = YAML::Load(c_str(false));
    scanKeys();
    unrecognizedKeywordCheck();
    checkParameters();
}

// Extract the key of each node as a list of strings.
// NOLINTNEXTLINE(misc-no-recursion)
static auto extractYAMLKeys(const YAML::Node& node,
                            std::vector<std::vector<std::string>>& keys,
                            std::vector<std::string>& cur_key) -> void
{
    if (node.IsMap()) {
        for (YAML::const_iterator it { node.begin() }; it != node.end(); ++it) {
            cur_key.push_back(it->first.as<std::string>());
            extractYAMLKeys(node[it->first], keys, cur_key);
            cur_key.pop_back();
        }
    } else {
        // Anything that is not a map means that we reached the end of
        // the tree of a configuration parameter (node). Store the
        // current key in the list and continue the search for the
        // next key.
        keys.push_back(cur_key);
    }
}

auto Settings::unrecognizedKeywordCheck() const -> void
{
    std::vector<std::vector<std::string>> all_keys {};
    std::vector<std::string> cur_key {};
    extractYAMLKeys(YAML::Clone(config), all_keys, cur_key);
    // Test all keys from the input file against all valid keys.
    for (const auto& key : all_keys) {
        bool found { false };
        for (const auto& valid_key : all_valid_keys) {
            if (key == valid_key) {
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "Warning: unrecognized input parameter: "
                      << Setting<bool> { key, false, "" }.keyToStr() << '\n';
        }
    }
}

auto Settings::c_str(const bool verbose) -> const char*
{
    do_dump = true;
    yaml_emitter.SetBoolFormat(YAML::YesNoBool);
    yaml_emitter.SetNullFormat(YAML::LowerNull);
    yaml_emitter.verbose = verbose;
    // Everything will be stored in a global map.
    yaml_emitter << YAML::BeginMap;
    scanKeys();
    do_dump = false;
    return yaml_emitter.c_str();
}

// Set a YAML node to a default (reference) value or keep the current
// value. It merges two YAML nodes by overriding the reference node if
// the corresponding key is present in the input node.
static auto setDefaultOrKeep(const YAML::Node& ref,
                             const YAML::Node& input) -> YAML::Node
{
    if (!input.IsMap()) {
        return input;
    }
    YAML::Node result { YAML::NodeType::Map };
    for (const auto& node : ref) {
        if (node.first.IsScalar()) {
            const std::string& key { node.first.Scalar() };
            if (input[key]) {
                result[node.first] = setDefaultOrKeep(node.second, input[key]);
                continue;
            }
        }
        result[node.first] = node.second;
    }
    return result;
}

auto Settings::getConfig() const -> std::string
{
    auto reduced_config { setDefaultOrKeep(default_config,
                                           YAML::Clone(config)) };
    YAML::Emitter out {};
    out.SetBoolFormat(YAML::YesNoBool);
    out.SetNullFormat(YAML::LowerNull);
    out << reduced_config;
    return out.c_str();
}

} // namespace tango
