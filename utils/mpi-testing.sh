#!/usr/bin/env bash

module load gnu openmpi/3.1.1 netcdf-fortran/4.4.4  pnetcdf/4.6.1
module load plumed/gnu
module load openblas

folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DINTERNAL_PLUMED=off -DWITH_EVB=On -DBUILDER="Gitlab Slave" && make -j10 && ctest --output-on-failure -j 2 -E TEST2[89]
