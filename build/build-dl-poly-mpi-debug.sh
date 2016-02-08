#!/usr/bin/env bash

COMP=intel
COMP=gcc
[[ -z $1 ]] && JOBS=10 || JOBS=$1
if [[ $COMP == "intel" ]]; then 
  module load intel/2016β1
  FC=ifort
  FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv"
else
  module load gnu-openmpi/1.8.3
  FC=gcc
  FFLAGS="-cpp -g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42"
fi

BUILD_DIR=/home/alin/playground/dl-poly-stfc-plumed/build/build-mpi-debug-$COMP
SOURCE_DIR=/home/alin/playground/dl-poly-stfc-plumed

[[ -d $BUILD_DIR ]] && rm -rf $BUILD_DIR 
mkdir -p $BUILD_DIR
pushd $BUILD_DIR
  cmake $SOURCE_DIR  -DWITH_MPI=y 
  make -j$JOBS
popd 
