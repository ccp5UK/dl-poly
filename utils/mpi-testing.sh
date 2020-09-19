#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6 netcdf-fortran/4.5.2  netcdf/4.7.3
module load plumed/gnu 
module load scalapack/2.1.0 dftb+/gcc/19.2
module load openblas/0.3.10

folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DINTERNAL_PLUMED=off -DWITH_EVB=On  && make -j10 && ctest --output-on-failure -j 2 -E TEST2[89]
