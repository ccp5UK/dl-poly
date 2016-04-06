#!/usr/bin/env bash
module load openmpi/gcc/1.10.1
module load netcdf/gcc/4.4.0
module load plumed/gcc/2.2.1
module load kim/gcc/1.7.2

mkdir build-mpi-all
pushd build-mpi-all
cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DWITH_NETCDF=ON 
make -j10 

