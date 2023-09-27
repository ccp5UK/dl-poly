#!/usr/bin/env bash

module load intel/2022a Python

mpr=`which mpirun`

rm -rf build-intel-testing
mkdir build-intel-testing
pushd build-intel-testing
FC=ifx cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWITH_ASSERT=ON -DMPIEXEC=$mpr -DWITH_EVB=ON   &&  make -j10 && ctest --output-on-failure -j 1 -E TEST2[89]
