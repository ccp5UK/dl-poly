#!/usr/bin/env bash

module purge
module load msmpi/7.1
rm -rf build-win32-msmpi
mkdir -p build-win32-msmpi
pushd build-win32-msmpi
FFLAGS="-O3 -static -fno-underscoring" cmake ../ -DCMAKE_TOOLCHAIN_FILE=../cmake/win32helper.cmake -DWITH_MPI=ON \
-DMPI_Fortran_INCLUDE_PATH="$MSMPIROOT/x86" -DMPI_Fortran_LIBRARIES="$MSMPIROOT/x86/libmsmpi.a"
make -j10
