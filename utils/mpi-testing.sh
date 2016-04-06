#!/usr/bin/env bash

module load openmpi/gcc/1.10.1
module load netcdf/gcc/4.4.0
module load plumed/gcc/2.2.1

mkdir build-mpi-testing
pushd build-mpi-testing
FFLAGS="-O3 -mtune=native" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DBUILDER="Gitlab Slave"
make -j10 
make test  


