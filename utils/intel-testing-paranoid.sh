#!/usr/bin/env bash

module load intel/2018.u3
mpr=`which mpirun` 
rm -rf build-intel-testing-p
mkdir build-intel-testing-p
pushd build-intel-testing-p
FC=ifort FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -qopt-report=5 -init=snan -init=arrays" cmake ../ -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave" -DWITH_PLUMED=off
make -j10
ctest -j 1 -E TEST28


