#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6

rm -rf build-mpi-kim
mkdir build-mpi-kim
pushd build-mpi-kim
rm -rf $HOME/101
cmake ../ -DWITH_KIM=ON -DCMAKE_INSTALL_PREFIX=$HOME/101
make -j10


