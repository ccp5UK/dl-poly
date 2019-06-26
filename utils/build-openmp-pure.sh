#!/usr/bin/env bash
module load gnu openblas
mkdir build-openmp-pure
pushd build-openmp-pure
cmake ../ -DWITH_MPI=OFF -DWITH_OPENMP=ON -DWITH_EVB=ON
make -j


