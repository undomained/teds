// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Extensions of the YAML library that are necessary to work with
// non-standard types of this project.

#pragma once

#include "constants.h"
#include "setting.h"
#include <yaml-cpp/yaml.h>
#include <string>
#include <unordered_map>
#include <spdlog/spdlog.h>
// Instruct YAML how to read values into non-standard types, e.g., how
// to read into a parameter of type std::optional.
namespace YAML {


// Structure to hold map data
struct MapResult {
    std::map<std::string, std::string> map;
    // Function to get the value by key
    // You can enter the key as: "nest1.nest2.key"
    std::string getVal(const std::string& key) const;
};

// Function to get all keys and values with nested keys as full path
// You can get the values in the map by using getVal() method
MapResult getMap(const YAML::Node& node, const std::string& prefix = "");


template <typename T>
struct convert<std::optional<T>>
{
    static auto decode(const Node& node, std::optional<T>& rhs) -> bool
    {
        // Null values are okay. Then the optional parameter remains unset.
        if (!node.IsNull()) {
            rhs.emplace(node.as<T>());
        }
        return true;
    }
};

template <>
struct convert<tango::ProcLevel>
{
    static auto decode(const Node& node, tango::ProcLevel& rhs) -> bool;
};

template <>
struct convert<tango::Unbin>
{
    static auto decode(const Node& node, tango::Unbin& rhs) -> bool;
};

} // namespace YAML

namespace tango {

// Instruct YAML how to print the value field of a parameter of type
// std::optional.
template <typename T>
auto operator<<(YAML::Emitter& out,
                const std::optional<T> value) -> YAML::Emitter&
{
    if (value) {
        out << value.value();
    } else {
        out << YAML::Null;
    }
    return out;
}

auto operator<<(YAML::Emitter& out,
                const ProcLevel proc_level) -> YAML::Emitter&;

auto operator<<(YAML::Emitter& out, const Unbin unbin) -> YAML::Emitter&;

// Extended Emitter to have a verbosity switch
class Emitter : public YAML::Emitter
{
public:
    bool verbose {};
};

// Instruct YAML how to print the verbose and non-verbose contents of
// a configuration parameter.
template <typename T>
static auto operator<<(Emitter& out, const Setting<T> setting) -> Emitter&
{
    out << YAML::Key << setting.yaml_keys.back();
    if (out.verbose) {
        out << YAML::Value;
        out << YAML::BeginMap;
        // NOLINTNEXTLINE(cppcoreguidelines-slicing)
        out << YAML::Key << "default" << YAML::Value << static_cast<T>(setting);
        out << YAML::Key << "type" << YAML::Value << setting.type;
        out << YAML::Key << "info" << YAML::Value << YAML::Literal
            << setting.info;
        out << YAML::EndMap;
    } else {
        // NOLINTNEXTLINE(cppcoreguidelines-slicing)
        out << YAML::Value << static_cast<T>(setting);
    }
    return out;
}

} // namespace tango
