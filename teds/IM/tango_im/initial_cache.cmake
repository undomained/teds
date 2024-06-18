# GNU compilers
set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")

# Flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -Wcast-align -Wformat -Winvalid-pch -Wmissing-declarations -Wmissing-include-dirs -Wconversion -Wredundant-decls -Wswitch-default -Wswitch-enum -pedantic" CACHE STRING "")

# Whether to build in Release or Debug configuration.
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# Where to install the CKD instrument model
set(CMAKE_INSTALL_PREFIX /usr/local/ckdmodel CACHE STRING "")

# Additional directories containing header files
set(INCLUDE_PATHS /usr/include CACHE STRING "")

# Library directories
set(LIBRARY_PATHS /usr/lib64 CACHE STRING "")

# Set Libraries
set(LIBRARIES netcdf_c++4 netcdf lapack CACHE STRING "")

# The IM additionally depends on the L1A-L1B processor
set(TANGO_L1B_PATH ../../L1AL1B/tango_l1b/build/ CACHE STRING "")

# Like with the L1A-L1B processor, point to where the dependencies are located.
set(SPDLOG_PATH /usr/local/spdlog CACHE STRING "")
set(YAML_CPP_PATH /usr/local/yaml-cpp CACHE STRING "")
set(POCKETFFT_PATH /usr/local/pocketfft CACHE STRING "")
