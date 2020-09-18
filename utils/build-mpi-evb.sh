#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6
module load openblas/0.3.10

mkdir build-mpi-pure-evb
pushd build-mpi-pure-evb
cmake ../ -DWITH_EVB=ON
make -j10


