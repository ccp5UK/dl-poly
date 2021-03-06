#!/usr/bin/env bash

module load gnu/7 openmpi/3.1.6
module load openblas/0.3.10
#module load plumed/gnu netcdf-fortran/4.4.4  pnetcdf/4.6.1

folder="build-mpi-testing-paranoid"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DWITH_PLUMED=Off -DWITH_EVB=On && make -j10 && ctest --output-on-failure -E TEST2[89]  --timeout 300
