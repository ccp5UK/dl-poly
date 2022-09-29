#!/usr/bin/env bash

module load intel/2022a

mkdir -p build-mpi-intel
pushd build-mpi-intel
FC=ifort cmake ../ -DMPI_Fortran_COMPILER=mpiifort 
make -j10

