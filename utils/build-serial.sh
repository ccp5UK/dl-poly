#!/usr/bin/env bash

mkdir build-serial
pushd build-serial
cmake ../ -DWITH_MPI=OFF 
make -j10 
