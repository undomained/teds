# Default variables for the release build configuration

set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG" CACHE STRING "")
set(LIBRARIES netcdf_c++4 netcdf hdf5_hl hdf5 z dl blas lapack CACHE STRING "")
