#!/usr/bin/env bash

module load gnu openmpi/3.0.0
module load kim/gcc

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON
make -j10


