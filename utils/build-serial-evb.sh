#!/usr/bin/env bash
module load gnu/7 openblas/0.3.10
mkdir build-serial
pushd build-serial
cmake ../ -DWITH_MPI=OFF -DWITH_EVB=On  
make -j10
