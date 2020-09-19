#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 scalapack/2.1.0 dftb+/gcc/19.2
build="build-mpi-dftbp"
rm -rf $build
mkdir $build
pushd $build
cmake ../ -DWITH_ASSERT=Off -DWITH_DFTBP=On -DCMAKE_BUILD_TYPE=Release
make -j10


