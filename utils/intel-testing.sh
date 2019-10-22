#!/usr/bin/env bash

module load intel/2018.u3
mpr=`which mpirun` 
rm -rf build-intel-testing
mkdir build-intel-testing
pushd build-intel-testing
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave"  && make -j10 && ctest -j 1 -E TEST28


