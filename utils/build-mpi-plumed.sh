#!/usr/bin/env bash

module load openmpi/gcc/1.10.1
module load netcdf/gcc/4.4.0
module load plumed/gcc/2.2.1

mkdir build-mpi-plumed
pushd build-mpi-plumed
cmake ../ -DWITH_PLUMED=ON 
make -j10 


