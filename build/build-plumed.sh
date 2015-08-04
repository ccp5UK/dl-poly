#!/usr/bin/env bash


module load gnu-openmpi/1.8.3
module list
CC=gcc
CXX=g++
PLUMED_VERSION=2.2b
PLUMED_SRC=plumed-${PLUMED_VERSION}.tgz
BUILD_FOLDER=build-$CC
PLUMED_INSTALL=/opt/plumed/$CC/$PLUMED_VERSION
PLUMED_MODULE=/opt/modules/plumed/$CC/$PLUMED_VERSION
CXXFLAGS="-O3 -mtune=native" 
CFLAGS="-O3 -mtune=native"
VMD_PLUGIN=/usr/lib64/vmd/1.9.2beta1.1427136773/plugins
VMD_INC="-I$VMD_PLUGIN/include -I$VMD_PLUGIN/LINUXAMD64/molfile "
LIBS="-lnetcdf -ltcl8.6 "
rm -rf $BUILD_FOLDER
mkdir -p $BUILD_FOLDER
pushd $BUILD_FOLDER
tar -xvf ../$PLUMED_SRC
pushd plumed-${PLUMED_VERSION}
patch -p1 < ../../0001-fix-template-issue-for-optimalAlignment.patch

./configure LDFLAGS="-L$VMD_PLUGIN/LINUXAMD64/molfile" CPPFLAGS="$VMD_INC" LIBS="$LIBS"  CXXFLAGS="$CXXFLAGS" CFLAGS="$CFLAGS" \
 --prefix=$PLUMED_INSTALL \
  --enable-external-molfile-plugins --enable-external-lapack \
 --enable-mpi --enable-openmp
make -j10
make install
popd 
popd 
mkdir -p $(dirname $PLUMED_MODULE)
cp $PLUMED_INSTALL/lib/plumed/src/lib/modulefile $PLUMED_MODULE
