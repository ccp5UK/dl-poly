#!/usr/bin/env bash
set -x
#module list
#COMPILER=intel
COMPILER=gnu
libsuffix=64
if [ $COMPILER == "gnu" ]; then
  module load openmpi/gcc
  export CC=gcc
  export CXX=mpic++
  export FC=gfortran
  export CXXFLAGS="-O3 -mtune=native" 
  export CFLAGS="-O3 -mtune=native"
else
  module load intel/2016
  export CC=icc
  export CXX=mpiicpc
  export FC=ifort
  export CXXFLAGS="-O3 -xHost -g" 
  export CFLAGS="-O3 -xHost -g"
fi
PLUMED_VERSION=2.2.1
BASE=/opt
PLUMED_FOLDER=plumed-${PLUMED_VERSION}
PLUMED_SRC=${PLUMED_FOLDER}.tgz
BUILD_FOLDER=build-$COMPILER
PLUMED_INSTALL=$BASE/plumed/$COMPILER/$PLUMED_VERSION
PLUMED_MODULE=$BASE/modules/plumed/$COMPILER/$PLUMED_VERSION
#VMD_PLUGIN=/usr/lib64/vmd/1.9.2beta1.1427136773/plugins
VMD_PLUGIN=
#VMD_INC="-I$VMD_PLUGIN/include -I$VMD_PLUGIN/LINUXAMD64/molfile "
VMD_INC=
LIBS="-lnetcdf -ltcl8.6 "
LDLF="-L$VMD_PLUGIN/LINUXAMD64/molfile"
rm -rf $BUILD_FOLDER
mkdir -p $BUILD_FOLDER
pushd $BUILD_FOLDER
tar -xf ../$PLUMED_SRC
pushd ${PLUMED_FOLDER}
#patch -p1 < ../../0001-fix-template-issue-for-optimalAlignment.patch
./configure LDFLAGS="$LDFL" CPPFLAGS="$VMD_INC" LIBS="$LIBS"  CXXFLAGS="$CXXFLAGS" CFLAGS="$CFLAGS" \
 --prefix=$PLUMED_INSTALL \
  --enable-external-molfile-plugins --enable-external-lapack \
 --enable-mpi --enable-openmp --enable-modules=crystallization
make -j10
make install
popd 
popd 
mkdir -p $(dirname $PLUMED_MODULE)
cp $PLUMED_INSTALL/lib${libsuffix}/plumed/modulefile $PLUMED_MODULE
