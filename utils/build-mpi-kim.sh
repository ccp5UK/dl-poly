#!/usr/bin/env bash

module load gnu openmpi/3.1.1
module load kim/gcc/2.1.3

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON -DINTERNAL_KIM=Off
make -j10


