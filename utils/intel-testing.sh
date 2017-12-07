#!/usr/bin/env bash

module load intel plumed/intel 
mpr=`which mpirun` 
mkdir build-intel-testing
pushd build-intel-testing
FC=ifort FFLAGS="-O1 " cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave" -DWITH_PLUMED=on
make -j10
make test


