#!/usr/bin/env bash

module load foss/2022a

mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../ 
make -j10


