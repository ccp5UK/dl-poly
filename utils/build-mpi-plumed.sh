#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 netcdf-fortran/4.5.2  netcdf/4.7.3
module load plumed/gnu
module load openblas

mkdir build-mpi-plumed
pushd build-mpi-plumed
cmake ../ -DWITH_PLUMED=ON -DINTERNAL_PLUMED=off -DWITH_EVB=on
make -j10


