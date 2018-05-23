#!/usr/bin/env bash

module load gnu/7 openmpi/3.0.0
mkdir build-mpi-pure
pushd build-mpi-pure
cmake ../
make -j10


