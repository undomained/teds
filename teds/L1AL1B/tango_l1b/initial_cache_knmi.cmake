# Initial cache variables for the SPEXone processor

# The syntax used here is
#  set(<variable> <value> CACHE <type> <docstring>) .
# You only need to edit the <value> parameter.

# GNU compilers
set(CMAKE_CXX_COMPILER /usr/lib64/ccache/g++ CACHE STRING "")
set(CMAKE_C_COMPILER /usr/lib64/ccache/gcc CACHE STRING "")

# GNU flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -Wcast-align -Wformat -Winvalid-pch -Wmissing-declarations -Wmissing-include-dirs -Wconversion -Wredundant-decls -Wswitch-default -Wswitch-enum -pedantic" CACHE STRING "")

# Intel compilers
#set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
#set(CMAKE_C_COMPILER icc CACHE STRING "")

# Intel flags
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -no-ipo -DNDEBUG" CACHE STRING "")
#set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g -debug all -traceback" CACHE STRING "")

# Whether to build in Release or Debug configuration.
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# Additional directories containing header files
set(INCLUDE_PATHS /usr/include CACHE STRING "")

# Linear algebra library directories
set(LIBRARY_PATHS /usr/lib64 CACHE STRING "")

# Linear algebra libraries
set(LIBRARIES libnetcdf_c++4.so libnetcdf.so libz.so libdl.so.2 libfftw3.so liblapack.so.3.11.0 CACHE STRING "")

# With MKL
#set(LIBRARIES netcdf_c++4 netcdf z dl fftw3 mkl_gf_lp64 mkl_sequential mkl_core CACHE STRING "")

# Whether to build with MPI support (usually automatically detected)
#set(USE_MPI ON CACHE BOOL "")

# Where to install the SPEXone processor
set(CMAKE_INSTALL_PREFIX /usr/local/ CACHE STRING "")