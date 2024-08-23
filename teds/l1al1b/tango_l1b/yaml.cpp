// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "yaml.h"

// Mapping between a process name and its corresponding entry in
// ProcLevel
const std::map<std::string, tango::ProcLevel> proc_level_to_enum {
    { "l1a", tango::ProcLevel::l1a },
    { "raw", tango::ProcLevel::raw },
    { "dark", tango::ProcLevel::dark },
    { "noise", tango::ProcLevel::noise },
    { "nonlin", tango::ProcLevel::nonlin },
    { "prnu", tango::ProcLevel::prnu },
    { "stray", tango::ProcLevel::stray },
    { "swath", tango::ProcLevel::swath },
    { "l1b", tango::ProcLevel::l1b },
    { "sgm", tango::ProcLevel::sgm },
};

// Mapping between a binning mode and its corresponding entry in Unbin
const std::map<std::string, tango::Unbin> unbin_to_enum {
    { "none", tango::Unbin::none },
    { "nearest", tango::Unbin::nearest },
    { "linear", tango::Unbin::linear },
    { "cubic", tango::Unbin::cubic },
};

namespace YAML {

auto convert<tango::ProcLevel>::decode(const Node& node,
                                       tango::ProcLevel& rhs) -> bool
{
    try {
        rhs = proc_level_to_enum.at(node.as<std::string>());
    } catch (const std::out_of_range&) {
        throw std::invalid_argument { "unknown process: "
                                      + node.as<std::string>() };
    }
    return true;
}

auto convert<tango::Unbin>::decode(const Node& node, tango::Unbin& rhs) -> bool
{
    try {
        rhs = unbin_to_enum.at(node.as<std::string>());
    } catch (const std::out_of_range&) {
        throw std::invalid_argument { "unknown binning mode: "
                                      + node.as<std::string>() };
    }
    return true;
}

// Implementation of the getVal function
std::string MapResult::getVal(const std::string& key) const {
    auto it = map.find(key);
    if (it != map.end()) {
        return it->second;
    } else {
        return ""; // Key not found
    }
}

// Helper function to recursively process YAML nodes
void processNode(const YAML::Node& node, const std::string& prefix, MapResult& result) {
    if (node.IsMap()) {
        for (const auto& item : node) {
            const std::string key = prefix.empty() ? item.first.as<std::string>() : prefix + "." + item.first.as<std::string>();
            const YAML::Node& value = item.second;

            // Recursively process nested nodes
            processNode(value, key, result);

            // Store key-value pairs
            if (value.IsScalar()) {
                const std::string valueStr = value.as<std::string>();
                result.map[key] = valueStr;
            }
        }
    } else if (node.IsSequence()) {
        for (std::size_t i = 0; i < node.size(); ++i) {
            const std::string key = prefix.empty() ? "[" + std::to_string(i) + "]" : prefix + "[" + std::to_string(i) + "]";
            const YAML::Node& item = node[i];

            // Recursively process sequence items
            processNode(item, key, result);

            // Store list items as values
            if (item.IsScalar()) {
                const std::string valueStr = item.as<std::string>();
                result.map[key] = valueStr;
            }
        }
    }
}

MapResult getMap(const YAML::Node& node, const std::string& prefix) {
    MapResult result;
    // If prefix is empty, process the entire node
    if (prefix.empty()) {
        processNode(node, "", result);
    } else {
        // Process only the subtree under the specified prefix
        std::string key = prefix;
        YAML::Node subNode = node;
        std::size_t pos = 0;
        std::string delimiter = ".";

        while ((pos = key.find(delimiter)) != std::string::npos) {
            std::string token = key.substr(0, pos);
            subNode = subNode[token];
            key.erase(0, pos + delimiter.length());
        }

        subNode = subNode[key]; // Process the final part
        processNode(subNode, prefix, result);
    }

    return result;

}



} // namespace YAML

namespace tango {

auto operator<<(YAML::Emitter& out,
                const ProcLevel proc_level) -> YAML::Emitter&
{
    const auto it { std::ranges::find_if(
      proc_level_to_enum,
      [proc_level](const auto& item) { return item.second == proc_level; }) };
    out << it->first;
    return out;
}

auto operator<<(YAML::Emitter& out, const Unbin unbin) -> YAML::Emitter&
{
    const auto it { std::ranges::find_if(
      unbin_to_enum,
      [unbin](const auto& item) { return item.second == unbin; }) };
    out << it->first;
    return out;
}

} // namespace tango
