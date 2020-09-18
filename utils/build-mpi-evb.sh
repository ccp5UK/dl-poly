#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6
module load openblas/0.3.10

mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../ -DWITH_EVB=ON
make -j10


