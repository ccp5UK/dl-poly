#!/usr/bin/env bash

module load openmpi/gcc/1.10.1
module load netcdf/gcc/4.4.0

mkdir build-mpi-netcdf
pushd build-mpi-netcdf
cmake ../ -DWITH_NETCDF=ON 
make -j10 


