#!/usr/bin/env bash

module load openmpi/gcc/1.10.1
module load kim/gcc/1.7.2

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON
make -j10


