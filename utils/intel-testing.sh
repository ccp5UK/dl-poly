#!/usr/bin/env bash

module load intel/2016.3 

mkdir build-intel-testing
pushd build-intel-testing
FC=ifort FFLAGS="-O3 -xHost" cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DBUILDER="Gitlab Slave"
make -j10
make test


