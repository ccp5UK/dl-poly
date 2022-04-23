#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 
module load plumed/gnu

mkdir build-mpi-plumed
pushd build-mpi-plumed
cmake ../ -DWITH_PLUMED=ON -DINTERNAL_PLUMED=off -DPLUMED_VERSION=2.4.2 
make -j10


