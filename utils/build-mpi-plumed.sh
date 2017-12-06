#!/usr/bin/env bash

module load openmpi/gcc
module load netcdf/gcc
module load plumed/gnu

mkdir build-mpi-plumed
pushd build-mpi-plumed
cmake ../ -DWITH_PLUMED=ON
make -j10


