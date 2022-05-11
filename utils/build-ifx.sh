#!/usr/bin/env bash

source /opt/intel/oneapi/setvars.sh

mkdir -p build-mpi-ifx
pushd build-mpi-ifx
FC=ifx cmake ../ -DMPI_BUILD_TYPE=Release
make -j10

