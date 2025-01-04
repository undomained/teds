// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Each user defined configuration parameter is stored in an instance
// of the Setting class. It is a templated class suitable for storing
// the values of primitive types as well as non-primitives such as
// vectors and strings. For a primitive type, the value is stored in
// the "value" field. For a non-primitive, the value can be accessed
// directly as if dealing with that type directly. This is achieved
// using multiple inheritance, e.g.
//
//   template <typename T>
//   class Setting<std::vector<T>>
//     : public SettingBase<T>
//     , public std::vector<T>
//
// Thus, for a vector, one can call things like size() directly on the
// Setting object. However, primitives cannot be used in multiple
// inheritance which is why their is stored separately in a member
// variable.

// The base class (SettingBase) stores the info and type fields which
// are not editable by the user but store helpful information. The
// yaml_keys field is a chain of node names (keys) used to find the
// parameter from a YAML configuration file. For instance, if the
// configuration file contains a section
//
//   prnu:
//     enabled: false
//
// then yaml_keys = { "prnu", "enabled" } for the "enabled"
// parameter. The last item in yaml_keys is the name of the parameter
// stored.
//
// For primitive types, the Setting class overloads several operators
// so that can be used without the need to reference the value field,
// e.g.
//
//   if (settings.prnu.enabled) {
//       ...
//
// In other words, we can use Setting<bool> as if it were bool. This
// is for convenience. The following has the same effect:
//
//   if (settings.prnu.enabled.value) {
//       ...
//
// Settings with non-primitive types such as vectors do not have the
// value field and can be operated as if they were that type
// (e.g. std::vector) due to multiple inheritance as described above.

#pragma once

#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace tango {

// Base class for storing some meta information about the setting and
// most importantly the yaml_keys field. Without it the setting could
// not be located in a configuration file. The actual value of the
// setting is defined in a child class.
template <typename T>
class SettingBase
{
public:
    std::vector<std::string> yaml_keys {};
    const std::string info {};
    std::string type {};

    // The default constructor defines an empty (unused) setting
    SettingBase() = default;
    SettingBase(const bool is_list,
                const std::vector<std::string>& yaml_keys,
                const std::string& info)
      : yaml_keys { yaml_keys }, info { info }
    {
        // Define a string representation of the value type. If the
        // setting is of the type list (e.g. std::vector), append
        // "list" to the type name.
        const std::string prefix { is_list ? " list" : "" };
        if constexpr (std::is_same_v<T, bool>) {
            type = "boolean" + prefix;
        } else if constexpr (std::is_same_v<T, int>) {
            type = "integer" + prefix;
        } else if constexpr (std::is_same_v<T, size_t>) {
            type = "unsigned integer" + prefix;
        } else if constexpr (std::is_same_v<T, double>) {
            type = "double (float64)" + prefix;
        } else if constexpr (std::is_same_v<T, std::string>) {
            type = "string" + prefix;
        } else if constexpr (std::is_enum_v<T>) {
            type = "string" + prefix;
        } else {
            throw std::domain_error {
                "type not supported by the Setting class. This can be fixed "
                "by introducing 'type' for this type in Setting."
            };
        }
    }
    // Convert a list of YAML keys into a string [a][b]...
    [[nodiscard]] auto keyToStr() const -> std::string
    {
        std::stringstream s {};
        for (const auto& key : yaml_keys) {
            s << '[' << key << ']';
        }
        return s.str();
    }
    ~SettingBase() = default;
};

// Setting class for primitive types
template <typename T>
class Setting : public SettingBase<T>
{
public:
    T value {};

    Setting() = default;
    Setting(const std::vector<std::string>& yaml_keys,
            const T value,
            const std::string& info)
      : SettingBase<T> { false, yaml_keys, info }, value { value }
    {}

    // Shortcuts that allow to bypass referencing value directly
    operator T() const { return value; }
    auto operator=(const T& value) -> Setting<T>&;
    auto operator*=(const T& value) -> Setting<T>&;
    auto operator/=(const T& value) -> Setting<T>&;

    ~Setting() = default;
};

template <typename T>
auto Setting<T>::operator=(const T& value) -> Setting<T>&
{
    this->value = value;
    return *this;
}

template <typename T>
auto Setting<T>::operator*=(const T& value) -> Setting<T>&
{
    this->value *= value;
    return *this;
}

template <typename T>
auto Setting<T>::operator/=(const T& value) -> Setting<T>&
{
    this->value /= value;
    return *this;
}

// Setting class for storing a list of values of some configuration
// parameter
template <typename T>
class Setting<std::vector<T>>
  : public SettingBase<T>
  , public std::vector<T>
{
public:
    Setting() = default;
    Setting(const std::vector<std::string>& yaml_keys,
            const std::vector<T> value,
            const std::string& info)
      : SettingBase<T> { true, yaml_keys, info }, std::vector<T> { value }
    {}
    auto operator=(const std::vector<T>& value) -> Setting<std::vector<T>>&;
};

// Make the assignment operator with a vector act on the std::vector
// base of the instance.
template <typename T>
auto Setting<std::vector<T>>::operator=(const std::vector<T>& value)
  -> Setting<std::vector<T>>&
{
    std::vector<T>* base { this };
    *base = value;
    return *this;
}

// Setting class for storing a list of lists of some configuration
// parameter
template <typename T>
class Setting<std::vector<std::vector<T>>>
  : public SettingBase<T>
  , public std::vector<std::vector<T>>
{
public:
    Setting() = default;
    Setting(const std::vector<std::string>& yaml_keys,
            const std::vector<T> value,
            const std::string& info)
      : SettingBase<T> { true, yaml_keys, info }
      , std::vector<std::vector<T>> { value }
    {}
    auto operator=(const std::vector<std::vector<T>>& value)
      -> Setting<std::vector<std::vector<T>>>&;
};

template <typename T>
auto Setting<std::vector<std::vector<T>>>::operator=(
  const std::vector<std::vector<T>>& value)
  -> Setting<std::vector<std::vector<T>>>&
{
    std::vector<std::vector<T>>* base { this };
    *base = value;
    return *this;
}

// Define a setting class for std::optional<T>. Here optional means
// that no reasonable default value could be defined in settings.h for
// that parameter but not necessarily that the user can omit defining
// that parameter in the configuration file. Similarly, there are
// parameters that are optional for the user but not defined as
// std::optinal<T> in settings.h because safe default values exist for
// those parameters.
//
// Note that the base class is called with just T which means that for
// std::optional<T> the string representation of type is still T as
// determined in the base class.
template <typename T>
class Setting<std::optional<T>>
  : public SettingBase<T>
  , public std::optional<T>
{
public:
    Setting() = default;
    Setting(const std::vector<std::string>& yaml_keys, const std::string& info)
      : SettingBase<T> { false, yaml_keys, info }, std::optional<T> {}
    {}
    auto operator=(const std::optional<T>& value) -> Setting<std::optional<T>>&;
    auto operator*=(const T& s) -> Setting<std::optional<T>>&;
    auto operator/=(const T& s) -> Setting<std::optional<T>>&;
};

// Same comment about the assignment operator as for std::vector above
template <typename T>
auto Setting<std::optional<T>>::operator=(const std::optional<T>& value)
  -> Setting<std::optional<T>>&
{
    std::optional<T>* base { this };
    *base = value;
    return *this;
}

template <typename T>
auto Setting<std::optional<T>>::operator*=(const T& value)
  -> Setting<std::optional<T>>&
{
    *this = this->value() * value;
    return *this;
}

template <typename T>
auto Setting<std::optional<T>>::operator/=(const T& value)
  -> Setting<std::optional<T>>&
{
    *this = this->value() / value;
    return *this;
}

// Define a setting class for strings
template <>
class Setting<std::string>
  : public SettingBase<std::string>
  , public std::string
{
public:
    Setting() = default;
    Setting(const std::vector<std::string>& yaml_keys,
            const std::string value,
            const std::string& info)
      : SettingBase<std::string> { false, yaml_keys, info }
      , std::string { value }
    {}
    // Same comment about the assignment operator as for std::vector above
    auto operator=(const std::string& value) -> Setting<std::string>&
    {
        std::string* base { this };
        *base = value;
        return *this;
    }
};

} // namespace tango
