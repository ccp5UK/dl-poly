#!/usr/bin/env bash

module load intel plumed/intel 
mpr=`which mpirun` 
rm -rf build-intel-testing-p
mkdir build-intel-testing-p
pushd build-intel-testing-p
FC=ifort FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -qopt-report=5 -init=snan -init=arrays" cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave" -DWITH_PLUMED=on
make -j10
make test


