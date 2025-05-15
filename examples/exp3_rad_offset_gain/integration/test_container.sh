#!/usr/bin/env sh

DIR=$(dirname $(realpath ${0}))

cd ${DIR} || exit 1

rm -f sample_l1b_offset.nc

docker run -ti \
  -v "$(pwd)":/data \
  -w /data \
  rad-offset-gain \
  --input sample_l1b.nc \
  --output sample_l1b_offset.nc \
  --offset 0.1
