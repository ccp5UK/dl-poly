#!/usr/bin/env bash

module load gnu openmpi/3.1.1 
#module load plumed/gnu netcdf-fortran/4.4.4  pnetcdf/4.6.1 

mkdir build-mpi-testing-paranoid
pushd build-mpi-testing-paranoid
FFLAGS="-g -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow" cmake ../ -DBUILD_TESTING=ON -DWITH_PLUMED=Off -DBUILDER="Gitlab Worker"
make -j10
ctest --timeout 200


