# Build instructions

Build instructions for the C++ code in the IM and L1AL1B blocks of the E2E simulator.

## Build IM program

```
cd teds/IM/tango_ckd_model
mkdir build && cd build
cmake -C ../initial_cache_knmi.cmake ..
make -j
```

This creates the executable `ckdmodel`.


## Build L1AL1B program

```
cd teds/L1AL1B/tango_l1b
mkdir build && cd build
cmake -C ../initial_cache_knmi.cmake ..
make -j
```

This creates the executable `tango_l1b`.