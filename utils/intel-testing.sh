#!/usr/bin/env bash

module load intel/2018.u3
module load gnu openblas
mpr=`which mpirun`
rm -rf build-intel-testing
mkdir build-intel-testing
pushd build-intel-testing
cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON  -DMPI_Fortran_COMPILER=mpiifort -DMPIEXEC=$mpr -DBUILDER="Gitlab Slave" -DWITH_EVB=ON   &&  make -j10 && ctest --output-on-failure -j 1 -E TEST2[89]
