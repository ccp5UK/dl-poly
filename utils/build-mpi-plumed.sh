#!/usr/bin/env bash

module load gnu openmpi/3.1.1 netcdf-fortran/4.4.4  pnetcdf/4.6.1
module load plumed/gnu

mkdir build-mpi-plumed
pushd build-mpi-plumed
cmake ../ -DWITH_PLUMED=ON
make -j10


