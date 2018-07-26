#!/usr/bin/env bash

module load gnu openmpi/3.0.0 netcdf-fortran/4.4.4  pnetcdf/4.6.1
module load plumed/gnu

mkdir build-mpi-testing
pushd build-mpi-testing
FFLAGS="-O3 -march=native -mtune=native" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DBUILDER="Gitlab Slave"
make -j10
ctest -j 2 


