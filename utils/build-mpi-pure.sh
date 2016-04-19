#!/usr/bin/env bash

module load openmpi/gcc/1.10.1
mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../
make -j10


