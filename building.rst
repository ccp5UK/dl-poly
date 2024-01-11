Building notes
==============

-  these notes are for building with `cmake <https://cmake.org>`__
-  you can pass options to cmake using **-DOPTION=value**. For a
   complete list of options inspect
   `cmake/DLPOLYBuildOptions.cmake <cmake/DLPOLYBuildOptions.cmake>`__
-  cmake -L <path to CMakeLists.txt> will show you a list of all
   available options.
-  explicit compiler specification can be achieved by using environment
   variable **FC** (eg. using Intel ifort *FC=ifort*)
-  compiler flags can be altered via **FFLAGS**, (eg *FFLAGS=“-O3
   -xHost”*)
-  one also can use **cmake-gui** or **ccmake** to setup the build
   options
-  to change the install path use **-DCMAKE_INSTALL_PREFIX=**
   (*-DCMAKE_INSTALL_PREFIX=$HOME/101/DL_POLY*)
-  automatic testing can be done after **DL_POLY_4** is built, using
   **make test**
-  to see all the tests available use **ctest -N**
-  to run one specific test use **ctest -R**
-  for a list of all supported targets **make help**
-  TODO check it works on Windows…

Standard MPI
------------

.. code:: sh

   git clone https://gitlab.com/ccp5/dl-poly.git
   cmake -S dl-poly -Bbuild-mpi-pure -DCMAKE_BUILD_TYPE=Release
   cmake --build build-mpi-pure
   cmake --install build-mpi-pure

older version of cmake

.. code:: sh

   mkdir build-mpi-pure
   pushd build-mpi-pure
   cmake ../ -DCMAKE_BUILD_TYPE=Release
   make -j10
   make install

-  will use whatever default MPI is found

*Intel Compilers - Intel MPI*

.. code:: sh

   FC=ifort cmake -S dl-poly -Bbuild-mpi-pure -DCMAKE_BUILD_TYPE=Release -DMPI_Fortran_COMPILER=mpiifort
   cmake --build build-mpi-pure
   cmake --install build-mpi-pure

*Intel Compilers - Some default mpi library, other than Intel MPI*

.. code:: sh

   FC=ifort cmake -S dl-poly -Bbuild-mpi-pure -DCMAKE_BUILD_TYPE=Release
   cmake --build build-mpi-pure
   cmake --install build-mpi-pure

Serial
------

.. code:: sh

   cmake -S dl-poly -Bbuild-serial -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF
   cmake --build build-serial
   cmake --install build-serial

*Intel Compilers*

.. code:: sh

   FC=ifort cmake -S dl-poly -Bbuild-serial -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF
   cmake --build build-serial
   cmake --install build-serial

Optimisation flags
------------------

-  gfortran/cray/intel use as options, flags are documented in
   `cmake/flags.cmake <cmake/flags.cmake>`__

.. code:: sh

   -DCMAKE_BUILD_TYPE=Release

.. code:: sh

   FFLAGS="-O3 -mtune=native"

-  Intel

.. code:: sh

   FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15"

-  If you plan to run the binary on a different type of a machine than
   you build it, check the manual of your compiler for the flags
   matching the *running machine*

Debugging, or when things go unexpected
---------------------------------------

-  gfortran/cray/intel use as options, flags are documented in
   `cmake/flags.cmake <cmake/flags.cmake>`__

.. code:: sh

   -DCMAKE_BUILD_TYPE=Debug

-  other compilers

.. code:: sh

   FFLAGS="desired flags" cmake ../

Building with KIM support
-------------------------

::

   cmake -S dl-poly -Bbuild-with-kim -DCMAKE_BUILD_TYPE=Release -DWITH_KIM=ON -DCMAKE_INSTALL_PREFIX=mypath
   cmake --build build-with-kim
   cmake --install build-with-kim

Building with PLUMED support
----------------------------

.. code:: sh

   cmake -S dl-poly -Bbuild-with-plumed -DCMAKE_BUILD_TYPE=Release -DWITH_PLUMED=ON -DCMAKE_INSTALL_PREFIX=mypath
   cmake --build build-with-plumed
   cmake --install build-with-plumed

building with FORD API Documentation
------------------------------------

.. code:: sh

   cmake -S dl-poly -Bbuild-with-ford -DDOCS_FORD=On
   cmake --build build-with-ford -- ford

FAQ
===

On Ubuntu machines
------------------

It was noticed that for some mpi implementations the linking stage
fails. You will see a lot of errors claiming undefined references to
MPI **solution**

.. code:: sh

   FC=mpif90 FFLAGS="-O3" cmake ../

Intel MPI
---------

Intel MPI Fortran wrapper breaks ifort preprocessing you will get an
error on the lines Len_trim(xxx) not supported or similar. **solution**
do not use FC=mpiifort
