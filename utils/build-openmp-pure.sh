#!/usr/bin/env bash
module load gnu/7
mkdir build-openmp-pure
pushd build-openmp-pure
cmake ../ -DWITH_MPI=OFF -DWITH_OPENMP=ON 
make -j


