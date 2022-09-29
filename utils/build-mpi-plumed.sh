#!/usr/bin/env bash

module load foss/2022a
module load PLUMED/2.8.0-foss-2022a

mkdir build-mpi-plumed
pushd build-mpi-plumed
FFLAGS="-fallow-argument-mismatch " cmake ../ -DWITH_PLUMED=ON -DINTERNAL_PLUMED=off -DPLUMED_VERSION=2.8.0 
make -j10


