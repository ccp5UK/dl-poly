#!/usr/bin/env bash

module load intel/2020a
mpr=`which mpirun`
rm -rf build-intel-testing-p
mkdir build-intel-testing-p
pushd build-intel-testing-p
FC=ifort cmake ../ -DBUILD_TESTING=ON  -DCMAKE_BUILD_TYPE=Debug -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DMPIEXEC_PREFLAGS="-check-mpi" -DWITH_PLUMED=off && make -j10 && ctest --output-on-failure -j 1 -E TEST2[89]
