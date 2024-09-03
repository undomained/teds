cd l1al1b/build
make -j
cd ../..
cd im/build
make -j
cd ../..
python run_E2E.py ../cfg/nitro/full_config.yaml im
