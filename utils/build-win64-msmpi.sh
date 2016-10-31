#!/usr/bin/env bash

module purge 
module load msmpi/7.1

rm -rf build-win64-msmpi
mkdir -p build-win64-msmpi
pushd build-win64-msmpi
FFLAGS="-O3 -static  -fno-underscoring" cmake ../ -DCMAKE_TOOLCHAIN_FILE=../cmake/win64helper.cmake -DWITH_MPI=ON \
-DMPI_Fortran_INCLUDE_PATH="$MSMPIROOT/x64" -DMPI_Fortran_LIBRARIES="$MSMPIROOT/x64/libmsmpi.a"
make -j10
