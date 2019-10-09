#!/usr/bin/env bash

module load gnu openmpi/3.1.1 netcdf-fortran/4.4.4  pnetcdf/4.6.1

mkdir build-mpi-plumed
pushd build-mpi-plumed
rm -rf $HOME/101
cmake ../ -DWITH_PLUMED=ON -DCMAKE_INSTALL_PREFIX=$HOME/101
make -j10


