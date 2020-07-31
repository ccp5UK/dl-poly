#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6
module load kim/gcc
module load openblas 

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON -DINTERNAL_KIM=Off -DWITH_EVB=On
make -j10


