#!/usr/bin/env bash 

COMPILER=intel
COMPILER=gnu
if [ $COMPILER == "gnu" ]; then
  module load openmpi/gcc/1.10.1 plumed/gnu/2.2.1
  export FC=gfortran 
  export FFLAGS="-cpp -O3"
  #export FFLAGS="-cpp -g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42"  cmake $SOURCE_DIR  -DWITH_MPI=y -DWITH_PLUMED=y
else
  module load intel/2016 plumed/intel/2.2.0
  export FC=ifort
  export FFLAGS="-fpp -O3 -g"
fi
VERSION=4.09β
DLPOLY_MODULE=/opt/modules/dlpoly/$COMPILER/$VERSION
BUILD_DIR=/home/alin/playground/dl-poly-alin-plumed/build/build-${COMPILER}-mpi-plumed-host
SOURCE_DIR=/home/alin/playground/dl-poly-alin-plumed

[[ -d $BUILD_DIR ]] && rm -rf $BUILD_DIR 
mkdir -p $BUILD_DIR
pushd $BUILD_DIR
  cmake $SOURCE_DIR  -DWITH_MPI=y -DWITH_PLUMED=y -DCMAKE_INSTALL_PREFIX=/opt/dlpoly/$COMPILER/$VERSION
  make -j10 VERBOSE=1
  make install 
  mkdir -p $(dirname $DLPOLY_MODULE)
  cp modulefile $DLPOLY_MODULE
popd 
