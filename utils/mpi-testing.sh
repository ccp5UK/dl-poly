#!/usr/bin/env bash

module load openmpi/gcc
module load netcdf/gcc
module load plumed/gnu

mkdir build-mpi-testing
pushd build-mpi-testing
FFLAGS="-O3 -march=native -mtune=native" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DBUILDER="Gitlab Slave"
make -j10
make test


