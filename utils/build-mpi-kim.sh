#!/usr/bin/env bash

module load openmpi/gcc
module load kim/gcc

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON
make -j10


