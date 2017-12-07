#!/usr/bin/env bash
set -x
#module list
#COMPILER=intel
COMPILER=gnu
libsuffix=64
if [ $COMPILER == "gnu" ]; then
  module load openmpi/gcc netcdf/4.5.0
  export CC=gcc
  export CXX=mpicxx
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
PLUMED_VERSION=2.3.3
[[ ! -f v$PLUMED_VERSION.tar.gz ]] && wget https://github.com/plumed/plumed2/archive/v$PLUMED_VERSION.tar.gz
BASE=/opt
PLUMED_FOLDER=plumed2-${PLUMED_VERSION}
PLUMED_SRC=v${PLUMED_VERSION}.tar.gz
BUILD_FOLDER=build-$COMPILER
PLUMED_INSTALL=$BASE/plumed/$COMPILER/$PLUMED_VERSION
PLUMED_MODULE=$BASE/modules/plumed/$COMPILER/$PLUMED_VERSION
#VMD_PLUGIN=/usr/lib64/vmd/1.9.2beta1.1427136773/plugins
VMD_PLUGIN=
#VMD_INC="-I$VMD_PLUGIN/include -I$VMD_PLUGIN/LINUXAMD64/molfile "
VMD_INC=
LIBS="-lnetcdf -ltcl8.6 "
#LDLF="-L$VMD_PLUGIN/LINUXAMD64/molfile"
rm -rf $BUILD_FOLDER
mkdir -p $BUILD_FOLDER
pushd $BUILD_FOLDER
tar -xf ../$PLUMED_SRC
pushd ${PLUMED_FOLDER}
#patch -p1 < ../../0001-fix-template-issue-for-optimalAlignment.patch
LDFLAGS="-L$NETCDF_DIR/lib64" LIBS="$LIBS"  CXXFLAGS="$CXXFLAGS" CFLAGS="$CFLAGS" ./configure \
 --prefix=$PLUMED_INSTALL \
  --enable-external-molfile-plugins --enable-external-lapack \
 --enable-mpi --enable-openmp --enable-modules=crystallization
make -j10
make install
popd 
popd 
mkdir -p $(dirname $PLUMED_MODULE)
cp $PLUMED_INSTALL/lib${libsuffix}/plumed/modulefile $PLUMED_MODULE
