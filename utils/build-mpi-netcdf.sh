#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 netcdf-fortran/4.5.2  netcdf/4.7.3 openblas
mkdir build-mpi-netcdf
pushd build-mpi-netcdf
cmake ../ -DWITH_NETCDF=ON -DWITH_EVB=ON
make -j10


