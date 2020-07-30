#!/usr/bin/env bash
module load gnu/7 openmpi/3.1.1

cp -r source source-mpi
pushd source-mpi
ln -s ../build/Makefile_MPI Makefile
make hpc


