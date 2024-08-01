#!/bin/bash

# Script to clone, build and install DISAMAR from GitLab

# usage:
# to install at specifed path: ./build_disamar.sh <path_to_install_location>
# to install in current folder: ./build_disamar.sh

# As default the program is compiled with gfortran

# Paper describing DISAMAR: https://doi.org/10.5194/gmd-15-7031-2022


loc="$1"
if [ -z "${loc}" ]
then
loc='.'
else
mkdir -p $loc
cd $loc
fi

git clone https://gitlab.com/KNMI-OSS/disamar/disamar.git/

cd disamar/src
pwd
exit

make
make install

cd ..
echo "DISAMAR installed in: $(pwd)/Disamar.exe"

cd "$(dirname "$0")"