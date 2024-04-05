#!/usr/bin/env bash
module load foss/2022a
mkdir build-serial-evb
pushd build-serial-evb
cmake ../ -DWITH_MPI=OFF -DWITH_EVB=On  
make -j10
