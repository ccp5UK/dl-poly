#!/usr/bin/env bash

module load foss/2022a Python
module load PLUMED/2.8.0-foss-2022a
module load kim-api/2.3.0-GCCcore-11.3.0

# Used in TEST186
kim-api-collections-management install user SNAP_ChenDengTran_2017_Mo__MO_698578166685_000

folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
export OMPI_MCA_rmaps_base_oversubscribe=true
FFLAGS="-fallow-argument-mismatch " cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWITH_PLUMED=ON -DPLUMED_VERSION=2.8.0 -DINTERNAL_PLUMED=off -DWITH_KIM=ON -DINTERNAL_KIM=Off -DWITH_EVB=On  && make -j10 && ctest --output-on-failure -j 2 -E TEST2[89]
