#!/usr/bin/env bash
module load openmpi/gcc/1.10.1

cp -r source source-mpi
pushd source-mpi
ln -s ../build/Makefile_MPI Makefile
make hpc


