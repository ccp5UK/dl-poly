#!/usr/bin/env bash

module load foss/2022a

mkdir build-mpi-pure-evb
pushd build-mpi-pure-evb
cmake ../ -DWITH_EVB=ON
make -j10


