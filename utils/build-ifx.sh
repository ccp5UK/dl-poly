#!/usr/bin/env bash

module load intel/2022a

mkdir -p build-mpi-ifx
pushd build-mpi-ifx
FC=ifx cmake ../ -DMPI_BUILD_TYPE=Release
make -j10

