to build using cmake
with gcc compilers with MPI and OpenMP support one needs to tell cmake -DMPI=on or y
-DOPENMP=y or on to activate netcdf support -DNETCDF=y

1. standard MPI only DLpoly

mkdir build
cd build
FFLAGS="-O3 -cpp" cmake ../ -DWITH_MPI=y

2. hybrid MPI and openmp DLPOLY

mkdir build
cd build
FFLAGS="-O3 -cpp" cmake ../ -DWITH_MPI=y -DWITH_OPENMP=y

3. serial DLPOLY
mkdir build
cd build
FFLAGS="-O3 -cpp" cmake ../ 

4. serial with OpenMP threads
mkdir build
cd build
FFLAGS="-O3 -cpp" cmake ../ -DWITH_OPENMP=y


intel compilers...
1. pure MPI
mkdir build
cd build
FC=mpiifort FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15" cmake ../ -DWITH_MPI=y

2. hybrid MPI and OpenMP
mkdir build
cd build
FC=mpiifort FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15" cmake ../ -DWITH_MPI=y -DWITH_OPENMP=y


For extra timing on node 0 add 
-DEXTRATIME=y to cmake
eg.
FC=mpiifort FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15" cmake ../ -DWITH_MPI=y -DWITH_OPENMP=y -DWITH_EXTRATIME=y

Xeon Phi:
for building native binaries for Xeon Phi pass to cmake -DWITH_PHI=y


Paranoid debug:
gfortran
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=42"
intel fortran
FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv"
known issues: 
1. FindMPI may trip with intel mpi
I_MPI_F90=ifort FC=ifort 
shall cure all issues.