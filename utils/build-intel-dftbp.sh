#!/usr/bin/env bash

module load intel/2019.0 dftb+/intel/19.2 
build="build-mpi-intel-dftb"
rm -rf $build
mkdir -p $build
pushd $build
FC=ifort cmake ../ -DMPI_Fortran_COMPILER=mpiifort -DSCALAPACK="MKL" -DCMAKE_BUILD_TYPE=Release -DWITH_ASSERT=Off -DWITH_DFTBP=On
make -j10 VERBOSE=1

