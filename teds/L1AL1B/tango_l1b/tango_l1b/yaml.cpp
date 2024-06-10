// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "yaml.h"

// Mapping between a process name and its corresponding entry in
// ProcLevel
const std::map<std::string, tango::ProcLevel> proc_level_str_to_enum {
    { "l1a", tango::ProcLevel::l1a },
    { "raw", tango::ProcLevel::raw },
    { "dark", tango::ProcLevel::dark },
    { "noise", tango::ProcLevel::noise },
    { "nonlin", tango::ProcLevel::nonlin },
    { "prnu", tango::ProcLevel::prnu },
    { "stray", tango::ProcLevel::stray },
    { "swath", tango::ProcLevel::swath },
    { "wave", tango::ProcLevel::wave },
    { "rad", tango::ProcLevel::rad },
    { "l1b", tango::ProcLevel::l1b },
};

namespace YAML {

auto convert<tango::ProcLevel>::decode(const Node& node, tango::ProcLevel& rhs)
  -> bool
{
    try {
        rhs = proc_level_str_to_enum.at(node.as<std::string>());
    } catch (const std::out_of_range&) {
        throw std::invalid_argument { "unknown process: "
                                      + node.as<std::string>() };
    }
    return true;
}

} // namespace YAML

namespace tango {

auto operator<<(YAML::Emitter& out, const ProcLevel proc_level)
  -> YAML::Emitter&
{
    const auto it { std::ranges::find_if(
      proc_level_str_to_enum,
      [proc_level](const auto& item) { return item.second == proc_level; }) };
    out << it->first;
    return out;
}

} // namespace tango
