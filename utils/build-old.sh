#!/usr/bin/env bash
module load gnu openmpi/3.0.0

cp -r source source-mpi
pushd source-mpi
ln -s ../build/Makefile_MPI Makefile
make hpc


