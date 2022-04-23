#!/usr/bin/env bash
module load gnu/7 openmpi/3.1.6 
module load plumed/gnu
module load kim/gcc
module load openblas/0.3.10

mkdir build-mpi-all
pushd build-mpi-all
cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DPLUMED_VERSION=2.4.2 -DINTERNAL_KIM=off -DINTERNAL_PLUMED=off -DWITH_EVB=on
make -j10

