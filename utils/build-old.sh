#!/usr/bin/env bash
module load openmpi/gcc

cp -r source source-mpi
pushd source-mpi
ln -s ../build/Makefile_MPI Makefile
make hpc


