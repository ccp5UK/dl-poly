#!/usr/bin/env bash

mkdir -p build-win64-serial
pushd build-win64-serial
FFLAGS="-O3 -static" cmake ../ -DCMAKE_TOOLCHAIN_FILE=../cmake/win64helper.cmake -DWITH_MPI=OFF -DWITH_EVB=ON
make -j10
