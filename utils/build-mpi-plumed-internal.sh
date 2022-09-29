#!/usr/bin/env bash

module load foss/2022a
mkdir build-mpi-plumed
pushd build-mpi-plumed
rm -rf $HOME/101
FFLAGS="-fallow-argument-mismatch " cmake ../ -DWITH_PLUMED=ON -DPLUMED_VERSION=2.8.0 -DCMAKE_INSTALL_PREFIX=$HOME/101
make -j10


