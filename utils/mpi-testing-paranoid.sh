#!/usr/bin/env bash

module load openmpi/gcc
module load netcdf/gcc
module load plumed/gnu

mkdir build-mpi-testing-paranoid
pushd build-mpi-testing-paranoid
FFLAGS="-g -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42 -ffpe-trap=invalid,zero,overflow -fdump-core" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DBUILDER="Gitlab Slave"
make -j10
make test


