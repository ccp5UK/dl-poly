#!/usr/bin/env bash

module load foss/2022a
mkdir build-serial
pushd build-serial
cmake ../ -DWITH_MPI=OFF  
make -j10
