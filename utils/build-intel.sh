#!/usr/bin/env bash

module load intel openblas
mkdir -p build-mpi-intel
pushd build-mpi-intel
FC=ifort cmake ../ -DMPI_Fortran_COMPILER=mpiifort -DWITH_EVB=ON
make -j10

