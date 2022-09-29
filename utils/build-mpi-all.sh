#!/usr/bin/env bash
module load foss/2022a 
module load kim-api/2.3.0-GCCcore-11.3.0
module load PLUMED/2.8.0-foss-2022a

mkdir build-mpi-all
pushd build-mpi-all
FFLAGS="-fallow-argument-mismatch " cmake ../ -DWITH_KIM=ON -DWITH_PLUMED=ON -DPLUMED_VERSION=2.8.0 -DINTERNAL_KIM=off -DINTERNAL_PLUMED=off -DWITH_EVB=on -DINTERACED=on
make -j10

