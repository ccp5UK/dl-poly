#!/usr/bin/env bash

module load gnu openmpi/3.1.1
#module load plumed/gnu netcdf-fortran/4.4.4  pnetcdf/4.6.1

folder="build-mpi-testing-paranoid"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DWITH_PLUMED=Off -DBUILDER="Gitlab Worker" && make -j10 && ctest --output-on-failure --timeout 200 -E TEST29
