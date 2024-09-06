cd l1al1b
rm -rf build
mkdir build
cd build
cmake -C ../initial_cache.cmake ..
make -j
cd ../..
cd im
rm -rf build
mkdir build
cd build
cmake -C ../initial_cache.cmake ..
make -j
cd ../..
python run_E2E.py ../cfg/nitro/full_config.yaml l1al1b
