#!/usr/bin/env bash

module load intel/2016.3 
mpr=`which mpirun` 
mkdir build-intel-testing
pushd build-intel-testing
FC=ifort FFLAGS="-O1 " cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave"
make -j10
make test


