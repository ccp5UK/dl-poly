#!/usr/bin/env bash

module load openmpi/gcc
module load netcdf/gcc

mkdir build-mpi-netcdf
pushd build-mpi-netcdf
cmake ../ -DWITH_NETCDF=ON
make -j10


