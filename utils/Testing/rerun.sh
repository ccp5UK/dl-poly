#!/usr/bin/env bash

module load openmpi/gcc 
DLROOT=$PWD/../../
pre=TEST
exe=$DLROOT/build-mpi/bin/DLPOLY.Z
mkdir -p new-data
pushd new-data
rm -rf *
for i in $(echo {01..25}); do
   test=$pre$i
   tar -xvf $DLROOT/data/$test.tar.gz
   pushd $test
   rm -f OUTPUT STATIS CFGMIN COLVAR OUTPUT.*
   mpirun -n 8 $exe
   rm -f REV*
   popd
done

