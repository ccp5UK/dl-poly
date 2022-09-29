#!/usr/bin/env bash

module load intel/2022a
mkdir -p build-mpi-intel-evb
pushd build-mpi-intel-evb
FC=ifort cmake ../ -DMPI_Fortran_COMPILER=mpiifort -DWITH_EVB=ON
make -j10

