#!/usr/bin/env bash

mkdir -p build-win32-serial
pushd build-win32-serial
FFLAGS="-O3 -static" cmake ../ -DCMAKE_TOOLCHAIN_FILE=../cmake/win32helper.cmake -DWITH_MPI=OFF 
make -j10
