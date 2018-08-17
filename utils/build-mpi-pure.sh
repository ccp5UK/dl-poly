#!/usr/bin/env bash

module load gnu openmpi/3.1.1
mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../
make -j10


