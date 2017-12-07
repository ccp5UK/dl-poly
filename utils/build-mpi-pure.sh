#!/usr/bin/env bash

module load openmpi/gcc
mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../
make -j10


