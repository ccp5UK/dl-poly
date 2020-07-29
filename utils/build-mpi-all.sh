#!/usr/bin/env bash
module load gnu/7 openmpi/3.1.6 netcdf-fortran/4.5.2  netcdf/4.7.3
module load plumed/gnu
module load kim/gcc

mkdir build-mpi-all
pushd build-mpi-all
cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DWITH_NETCDF=ON -DINTERNAL_KIM=off -DINTERNAL_PLUMED=off
make -j10

