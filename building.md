# Building notes
* these notes are for building with [**cmake**](https://cmake.org)
* you can pass options to cmake using **-DOPTION=value**. For a complete list of options inspect [cmake/DLPOLYBuildOptions.cmake](cmake/DLPOLYBuildOptions.cmake)
* cmake -L <path to CMakeLists.txt> will show you a list of all available options.
* explicit compiler specification can be achieved by using environment variable **FC** (eg. using Intel ifort *FC=ifort*)
* compiler flags can be altered via **FFLAGS**, (eg *FFLAGS="-O3 -xHost"*)
* one also can use **cmake-gui** or **ccmake** to setup the build options
* to change the install path use **-DCMAKE_INSTALL_PREFIX=<path>** (*-DCMAKE_INSTALL_PREFIX=$HOME/101/DL_POLY*)
* automatic testing can be done after **DL_POLY_4** is built, using **make test**
* to see all the tests available use **ctest -N**
* to run one specific test use **ctest -R <TESTNAME>**
* for a list of all supported targets **make help**
* TODO check it works on Windows...

## Standard MPI
```sh
mkdir build-mpi-pure
pushd build-mpi-pure
FFLAGS="-O3" cmake ../
make -j10
make install
```
* will use whatever default MPI is found

*Intel Compilers - Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DMPI_Fortran_COMPILER=mpiifort 
```
*Intel Compilers - Some default mpi library, other than Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ 
```

## Hybrid MPI and OpenMP
```sh
mkdir build-mpi-openmp
pushd build-mpi-openmp
FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON
make -j10
make install
```
*Intel Compilers - Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DMPI_Fortran_COMPILER=mpiifort
```

## Serial
```sh
mkdir build-serial
pushd build-serial
FFLAGS="-O3" cmake ../ -DWITH_MPI=OFF
```
*Intel Compilers*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_MPI=OFF
```

## Serial with OpenMP threads
```sh
mkdir build-openmp
pushd build-openmp
FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DWITH_MPI=OFF
```
*Intel Compilers*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DWITH_MPI=OFF
```
## Xeon Phi
```sh
mkdir build-xeonphi
pushd build-xeonphi
FC=ifort FFLAGS="-O3 " cmake ../ -DWITH_PHI=ON -DWITH_MPI=ON
```

## Optimisation flags
* gfortran

```sh
FFLAGS="-O3 -mtune=native"
```

* Intel

```sh
FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15"
```

* If you plan to run the binary on a different type of a machine than you build it, check the manual of your compiler
for the flags matching the _running machine_

## Debugging, or when things go merdre
* gfortran

```sh
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=35 -ffpe-trap=invalid,zero,overflow -fdump-core"
```
* Intel

```sh
FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv"
```

## Building with NETCDF support
```sh
mkdir build-mpi-netcdf
pushd build-mpi-netcdf
FFLAGS="-O3" cmake ../ -DWITH_NETCDF=ON
make -j10
make install
```

## Building with KIM support
```
mkdir build-mpi-kim
pushd build-mpi-kim
FFLAGS="-O3" cmake ../ -DWITH_KIM=ON
make -j10
make install
```

## Building with PLUMED support
```sh
mkdir build-mpi-plumed
pushd build-mpi-plumed
FFLAGS="-O3" cmake ../ -DWITH_PLUMED=ON
make -j10
make install
```


# FAQ

## On Ubuntu machines

it was noticed that for some mpi implementations the linking stage fails. You will see a lot of errors claiming undefined references to MPI_*
**solution**

```sh
FC=mpif90 FFLAGS="-O3" cmake ../
```

## Intel MPI

Intel MPI Fortran wrapper breaks ifort preprocessing
you will get an erro on the lines Len_trim(xxx) not supported or similar.
**solution**
do not use FC=mpiifort
