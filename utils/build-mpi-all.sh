#!/usr/bin/env bash
module load openmpi/gcc
module load netcdf/gcc
module load plumed/gnu
module load kim/gcc

mkdir build-mpi-all
pushd build-mpi-all
cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DWITH_NETCDF=ON
make -j10

