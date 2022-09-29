#!/usr/bin/env bash

module load foss/2022a kim-api/2.3.0-GCCcore-11.3.0 

mkdir build-mpi-kim
pushd build-mpi-kim
cmake ../ -DWITH_KIM=ON -DINTERNAL_KIM=Off 
make -j10


