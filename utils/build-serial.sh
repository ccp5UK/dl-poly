#!/usr/bin/env bash
module load gnu openblas
mkdir build-serial
pushd build-serial
cmake ../ -DWITH_MPI=OFF  -DWITH_EVB=ON
make -j10
