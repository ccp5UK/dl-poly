#!/usr/bin/env bash
module load foss/2022a
mkdir build-openmp-pure
pushd build-openmp-pure
cmake ../ -DWITH_MPI=OFF -DWITH_OPENMP=ON 
make -j


