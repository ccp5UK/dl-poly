#!/usr/bin/env bash

module load intel/2020a
mpr=`which mpirun`
rm -rf build-intel-testing
mkdir build-intel-testing
pushd build-intel-testing
FC=ifort cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr &&  make -j10 && ctest --output-on-failure -j 1 -E TEST2[89]
