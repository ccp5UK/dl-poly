#!/usr/bin/env bash
module load gnu openmpi/3.1.1 netcdf-fortran/4.4.4  pnetcdf/4.6.1
module load plumed/gnu
module load kim/gcc/2.1.3
module load openblas

mkdir build-mpi-all
pushd build-mpi-all
cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DWITH_NETCDF=ON -DINTERNAL_KIM=off -DINTERNAL_PLUMED=off -DWITH_EVB=on
make -j10

