#!/usr/bin/env bash

module load gnu/7 openmpi/3.0.0 
#module load plumed/gnu netcdf-fortran/4.4.4  pnetcdf/4.6.1 

mkdir build-mpi-testing-paranoid
pushd build-mpi-testing-paranoid
FFLAGS="-g -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42 -ffpe-trap=invalid,zero,overflow -fdump-core" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=Off -DBUILDER="Gitlab Slave"
make -j10
make test


