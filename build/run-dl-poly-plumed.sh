#!/usr/bin/env bash 

module load gnu-openmpi/1.8.3 plumed/gcc/2.2b

DL_POLY=/home/alin/playground/dl-poly-stfc-plumed/build/build-mpi-plumed-host/bin/DLPOLY.Z
TEST_DIR=/home/alin/playground/dl-poly-stfc-plumed/build/plumed-test

[[ -d $TEST_DIR ]] && rm -rf $TEST_DIR
mkdir -p $TEST_DIR
cp CONFIG $TEST_DIR/
cp FIELD $TEST_DIR/
cp CONTROL $TEST_DIR/
cp plumed-biased.dat $TEST_DIR/
cp plumed-metadyn.dat $TEST_DIR/
pushd $TEST_DIR
  mpirun -n 1 $DL_POLY 
popd
