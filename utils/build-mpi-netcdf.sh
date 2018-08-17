#!/usr/bin/env bash

module load gnu openmpi/3.1.1 netcdf-fortran/4.4.4  pnetcdf/4.6.1 
mkdir build-mpi-netcdf
pushd build-mpi-netcdf
cmake ../ -DWITH_NETCDF=ON
make -j10


