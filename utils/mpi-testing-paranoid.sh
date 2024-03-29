#!/usr/bin/env bash

module load foss/2022a Python


folder="build-mpi-testing-paranoid"
rm -rf $folder && mkdir $folder && pushd $folder
export OMPI_MCA_rmaps_base_oversubscribe=true
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -DWITH_PLUMED=Off -DWITH_EVB=On && make -j10 && ctest --output-on-failure -E TEST2[89]  --timeout 300
