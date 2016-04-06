#!/usr/bin/env bash

module load intel/2016.1
FC=ifort cmake ../ -DMPI_Fortran_COMPILER=mpiifort
make -j8

