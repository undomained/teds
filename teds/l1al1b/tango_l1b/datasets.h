// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include <set>
#include <unordered_map>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>
#include "ckd.h"
#include "binning_table.h"
#include "l1.h"

namespace tango {

class DataContainer {
private:
    std::string c_name;
    std::variant<L1,CKD,BinningTable> this_container;

public:
    /// Constructor.
    template<typename T>
    DataContainer(
            std::string container_name,
            T container
            ) : c_name{ std::move(container_name) }, this_container{container}{}

    std::string get_name(){
        return c_name;
    }
    template <typename T>
    T* get(){
        return &std::get<T>(this_container);
    }
    template <typename T>
    T const* get() const {
        return &std::get<T>(this_container);
    }

    friend bool operator <(const DataContainer& lhs, const DataContainer& rhs) { //sort using c_name 
      return lhs.c_name < rhs.c_name;
   }
    friend bool operator <(const std::string& name, const DataContainer& rhs) { //sort using c_name 
      return name < rhs.c_name;
   }
    friend bool operator <(const DataContainer& lhs, const std::string& name) { //sort using c_name 
      return lhs.c_name < name;
   }
};
class Dataset {
private:
    std::set<DataContainer,std::less<>> datasets;

public:
    /// Constructor.
    Dataset() = default;

    /// Destructor.
    ~Dataset() = default;

    template <typename T> 
    void add_container(const std::string c_name, T container){
        datasets.emplace(DataContainer{c_name, container});
    };

    template <typename T> 
    T const& get_container(const std::string c_name) const{
        auto i = datasets.find( c_name );
        if ( i == datasets.end() ) {
        	spdlog::error("Container {} does not exist in datasets!!!", c_name);
                throw std::runtime_error("rquested container does not exist");
        }
        T const* ptr = i->template get<T>();
        if (! ptr ) { spdlog::error("container of wrong type",c_name); 
           throw std::runtime_error("wrong type of container");
        }
        return *ptr;
    };

};



} // namespace tango
