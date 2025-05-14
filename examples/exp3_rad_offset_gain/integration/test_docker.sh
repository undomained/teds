#!/usr/bin/env sh

rm -f sample*nc

python3 nc_generator.py || exit 1

docker run -ti \
  -v "$(pwd)":/data \
  -w /data \
  rad-offset-gain \
  --input sample_l1b.nc \
  --output sample_l1b_offset.nc \
  --offset 0.1
