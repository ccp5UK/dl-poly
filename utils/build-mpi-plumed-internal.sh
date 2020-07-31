#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 netcdf-fortran/4.5.2  netcdf/4.7.3

mkdir build-mpi-plumed
pushd build-mpi-plumed
rm -rf $HOME/101
cmake ../ -DWITH_PLUMED=ON -DCMAKE_INSTALL_PREFIX=$HOME/101
make -j10


