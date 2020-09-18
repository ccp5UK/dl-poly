#!/usr/bin/env bash
module load gnu/7 openblas/0.3.10
mkdir build-serial-evb
pushd build-serial-evb
cmake ../ -DWITH_MPI=OFF -DWITH_EVB=On  
make -j10
