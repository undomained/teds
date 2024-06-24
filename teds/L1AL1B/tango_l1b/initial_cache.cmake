# Set compilers
set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")

# Set compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -Wcast-align -Wformat -Winvalid-pch -Wmissing-declarations -Wmissing-include-dirs -Wconversion -Wredundant-decls -Wswitch-default -Wswitch-enum -pedantic" CACHE STRING "")

# Whether to build in Release or Debug configuration.
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# Additional directories containing header files
set(INCLUDE_PATHS /usr/include CACHE STRING "")

# Linear algebra library directories
set(LIBRARY_PATHS /usr/lib64 CACHE STRING "")

# Linear algebra libraries
set(LIBRARIES netcdf_c++4 netcdf lapack CACHE STRING "")

# Where to install
set(CMAKE_INSTALL_PREFIX /usr/local/tango CACHE STRING "")

# If installed, point to where the dependencies are located. This is
# not necessary if they are in a standard location installed by your
# package manager.
set(SPDLOG_PATH /usr/local/spdlog CACHE STRING "")
set(YAML_CPP_PATH /usr/local/yaml-cpp CACHE STRING "")
set(POCKETFFT_PATH /usr/local/pocketfft CACHE STRING "")
