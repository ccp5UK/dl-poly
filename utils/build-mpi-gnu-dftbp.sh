#!/usr/bin/env bash

module load gnu openmpi/3.1.1 scalapack/2.0.2 dftb+/gcc/19.2
build="build-mpi-dftbp"
rm -rf $build
mkdir $build
pushd $build
cmake ../ -DWITH_ASSERT=Off -DWITH_DFTBP=On -DCMAKE_BUILD_TYPE=Release
make -j10


