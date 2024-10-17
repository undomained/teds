rm -rf build
mkdir build
cd build
cmake -C ../initial_cache.cmake ..
make -j
