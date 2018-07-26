#!/usr/bin/env bash

module load intel 
mpr=`which mpirun` 
mkdir build-intel-testing
pushd build-intel-testing
FC=ifort FFLAGS="-O3 -g -traceback" cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave" 
make -j10
ctest -j 2 -E TEST28


