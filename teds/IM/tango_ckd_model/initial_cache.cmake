# Initial cache variables for the SPEXone processor

# The syntax used here is
#  set(<variable> <value> CACHE <type> <docstring>) .
# You only need to edit the <value> parameter.

# GNU compilers
set(CMAKE_CXX_COMPILER /usr/lib64/ccache/g++ CACHE STRING "")
set(CMAKE_C_COMPILER /usr/lib64/ccache/gcc CACHE STRING "")

# Flags
set(CMAKE_CXX_FLAGS -std=c++17 CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -Wcast-align -Wformat -Winvalid-pch -Wmissing-declarations -Wmissing-include-dirs -Wconversion -Wredundant-decls -Wswitch-default -Wswitch-enum -pedantic" CACHE STRING "")

# Whether to build in Release or Debug configuration.
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# Whether to build with MPI support
#set(USE_MPI ON CACHE BOOL "")

# Where to install the CKD instrument model
set(CMAKE_INSTALL_PREFIX /usr/local/ckdmodel CACHE STRING "")

# Keyword for automatic detection of a linear algebra library. Allowed
# values are:
#   mkl
#   none*
# *choosing none explicitly forces to use the fallback option
#set(LINALG_FLAVOR mkl CACHE STRING "")

# Additional directories containing header files
set(INCLUDE_PATHS /usr/include  /net/pc230061/nobackup/users/leune/bin/netcdf-cxx4/include/ CACHE STRING "")

# Library directories
set(LIBRARY_PATHS /usr/lib64 /net/pc230061/nobackup/users/leune/bin/netcdf-cxx4/lib64 CACHE STRING "")

# Set Libraries
#set(LIBRARIES libnetcdf_c++.so.4.2.0 libnetcdf.so libz.so libdl.so.2 libfftw3.so liblapack.so.3.11.0 CACHE STRING "")
set(LIBRARIES libnetcdf-cxx4.so libnetcdf.so libz.so libdl.so.2 libfftw3.so liblapack.so.3.11.0 CACHE STRING "")

