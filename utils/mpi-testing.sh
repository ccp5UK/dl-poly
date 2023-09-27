#!/usr/bin/env bash

module load foss/2022a Python
module load PLUMED/2.8.0-foss-2022a


folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
FFLAGS="-fallow-argument-mismatch " cmake ../ -DCMAKE_BUILD_TYPE=Release -DWITH_ASSERT=ON -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DPLUMED_VERSION=2.8.0 -DINTERNAL_PLUMED=off -DWITH_EVB=On  && make -j10 && ctest --output-on-failure -j 2 -E TEST2[89]
