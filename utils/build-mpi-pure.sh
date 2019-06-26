#!/usr/bin/env bash

module load gnu openmpi/3.1.1 openblas
mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../ -DWITH_EVB=ON
make -j10


