#!/usr/bin/env bash

mkdir build-openmp-pure
pushd build-openmp-pure
cmake ../ -DWITH_MPI=OFF -DWITH_OPENMP=ON
make -j


