#!/usr/bin/env bash

module load intel/2018.u3
mpr=`which mpirun` 
rm -rf build-intel-testing-p
mkdir build-intel-testing-p
pushd build-intel-testing-p
FC=ifort cmake ../ -DBUILD_TESTING=ON  -DCMAKE_BUILD_TYPE=Debug -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DMPIEXEC_PREFLAGS="-check-mpi" -DBUILDER="Gitlab Slave" -DWITH_PLUMED=off && make -j10 && ctest -j 1 -E TEST2[89]


