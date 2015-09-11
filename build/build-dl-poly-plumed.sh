#!/usr/bin/env bash 

module load gnu-openmpi/1.8.3 plumed/gcc/2.2b

BUILD_DIR=/home/alin/playground/dl-poly-stfc-plumed/build/build-mpi-plumed-host
SOURCE_DIR=/home/alin/playground/dl-poly-stfc-plumed

[[ -d $BUILD_DIR ]] && rm -rf $BUILD_DIR 
mkdir -p $BUILD_DIR
pushd $BUILD_DIR
  FC=gfortran FFLAGS="-cpp -g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42"  cmake $SOURCE_DIR  -DWITH_MPI=y -DWITH_PLUMED=y
  make -j10
popd 
