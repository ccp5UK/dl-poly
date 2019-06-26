#!/usr/bin/env bash

module load gnu openmpi/3.1.1
module load kim/gcc/2.0.0
module load openblas

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON -DWITH_EVB=ON
make -j10


