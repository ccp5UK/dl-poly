DL_POLY_4.01 CUDA+OpenMP PORT - an ICHEC-STFC Collaboration
===========================================================

Developed by Christos Kartsaklis and Ruairi Nestor
(ruairi.nestor@ichec.ie) in collaboration with
Ilian Todorov (ilian.todorov@stfc.ac.uk).


PORT STATUS
===========
GPU assisted compuation is only available to the following subset
of the DL_POLY_4 features/functionality/routines:

* link_cell_pairs - building the Verlet Neighbour List using Link
                    Cells
* metal_ld_compute - provided metal interactions are non-EAM
                     (keypot!=0)
* metal_forces - provided the default/tabulated calculation is
                 assumed (i.e. no "metal direct" evaluation is
                 opted for in CONTROL, ld_met controlled)
* vdw_forces - provided the default/tabulated calculation is
               assumed (i.e. no "vdw direct" evaluation is opted
               for in CONTROL, ld_vdw controlled)
* ewald_real_forces - real space Ewald electrostatics
* ewald_spme_forces - reciprocal space Ewald electrostatics
                      without the 3D DFFT (DaFT/GPFA)
* constraints_shake_vv and constraints_shake_lfv
* images (in numeric_container) - only for image conditions 2 & 3

INSTALLATION INSTRUCTIONS
=========================
Copy the Makefile_CUDA from the CUDA sub-directory as a Makefile
into the DL_POLY_4 source directory above.  This Makefile may need
to be configured for the system on which the code will be compiled
so that libraries and include files can be found by the compilers.

COMPILATION INFORMATION
=======================
The code is known to compile with either of CUDA SDK 2.3/3.1 and for
compute capabilities 1.3 and 2.0.  If OpenMP is used, a GCC and
FORTRAN90 compiler that support OpenMP are required to build the
sources.

If single precision is desired, as soon as the FORTRAN90 sources'
precision has been modified in kinds_f90.f90, change
CFG_DOUBLE_PRECISION to 0 in dl_poly_cu.h and recompile.

Depending on the compute cability of the GPU (1.3 and 2.0 are
supported currently), the CFG_COMPUTE_{MAJOR, MINOR} constants have
to be set accordingly in the Makefile.

Concurrent execution of accelerated functions on the host and device
is enabled by default. To disable concurrent execution
(i.e. so that accelerated functions execute on the device only), the
preprocessor vairable CFG_OVERLAP_WITH_HOST should be set to 0 in
dl_poly_cu.h before compilation.

Finally, CUDA (& OpenMP, if present) acceleration for a particular
function can be disabled by setting the appropriate
dl_poly_cuda_offload_* function in dl_poly_init_cu.cu.  This can
be helpful during debugging.

EXECUTION INSTRUCTIONS & EXAMPLES
=================================
When the DLPOLY.Z.cu executable is invoked, an initialiser in
dl_poly_init_cu.cu is invoked form dl_poly.f90.  The initialiser
binds the MPI processes running on each node to the available
CUDA-enabled devices attached to that node in round-robin fashion.
Information regarding the MPI process affinity and device binding
is printed to the standard output.

Please NOTE that the affinity of the MPI processes is not decided
by the initializer - it is the responsibility of the user to
specify this when the executable is run.  The documentation for
the particular MPI implementation in use should be consulted in
order to achieve this.

Examples of the MPI process-device binding and process affinity
output (on a system with 8-core nodes, 2 CUDA devices per node,
and using MVAPICH2) are:

Example 1: 2 MPI processes, 1 node:
===
./cleanup; mpiexec -n 2 DLPOLY.Z.cu
....
Bound MPI process 0 (pid=14299; affined to CPU(s)  0; 1 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 1 (pid=14300; affined to CPU(s)  1; 1 OpenMP thread(s)) to device 1@stoney52

Example 2: 2 MPI process, 1 node, 4 threads each:
===
export OMP_NUM_THREADS=4
./cleanup; mpiexec -n 2 DLPOLY.Z.cu

Bound MPI process 0 (pid=14319; affined to CPU(s)  0; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 1 (pid=14320; affined to CPU(s)  1; 4 OpenMP thread(s)) to device 1@stoney52

Example 3: 8 MPI processes, 1 node (GPU over-subscription):
===
./cleanup ; mpiexec -n 8

Bound MPI process 0 (pid=14345; affined to CPU(s)  0; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 1 (pid=14346; affined to CPU(s)  1; 4 OpenMP thread(s)) to device 1@stoney52
Bound MPI process 2 (pid=14347; affined to CPU(s)  2; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 3 (pid=14348; affined to CPU(s)  3; 4 OpenMP thread(s)) to device 1@stoney52
Bound MPI process 4 (pid=14349; affined to CPU(s)  4; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 5 (pid=14350; affined to CPU(s)  5; 4 OpenMP thread(s)) to device 1@stoney52
Bound MPI process 6 (pid=14351; affined to CPU(s)  6; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 7 (pid=14352; affined to CPU(s)  7; 4 OpenMP thread(s)) to device 1@stoney52

** WARNING: The number of MPI processes (8) on node stoney52 is greater than the number of devices (2) on that node
   Ideally, the number of MPI processes on a given node should match the number of available devices on that node.


Example 4: 2 MPI processes, 1 node, 4 threads each, explicit affinity (set using the
MVAPICH2 environment variable, MV2_CPU_MAPPING):
===
export OMP_NUM_THREADS=4
export MV2_CPU_MAPPING="0,2,4,6:1,3,5,7"
./cleanup ; mpiexec -n 2 ./DLPOLY.Z.cu
….
Bound MPI process 0 (pid=14409; affined to CPU(s)  0 2 4 6; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 1 (pid=14410; affined to CPU(s)  1 3 5 7; 4 OpenMP thread(s)) to device 1@stoney52


Example 5: 4 MPI processes, 2 nodes, 4 threads each, explicit affinity (set using the
MVAPICH2 environment variable MV2_CPU_MAPPING):
===
export OMP_NUM_THREADS=4
export MV2_CPU_MAPPING="0,2:4,6:1,3:5,7"
./cleanup ; mpiexec -n 4 -npernode 2 ./DLPOLY.Z.cu
….
Bound MPI process 0 (pid=14627; affined to CPU(s)  0 2 4 6; 4 OpenMP thread(s)) to device 0@stoney52
Bound MPI process 1 (pid=14628; affined to CPU(s)  1 3 5 7; 4 OpenMP thread(s)) to device 1@stoney52
Bound MPI process 2 (pid=8581; affined to CPU(s)  0 2 4 6; 4 OpenMP thread(s)) to device 0@stoney51
Bound MPI process 3 (pid=8582; affined to CPU(s)  1 3 5 7; 4 OpenMP thread(s)) to device 1@stoney51

FURTHER NOTES
=============
Apart from the resource organisation, a number of shell variables
modify the runtime behavior; these are documented inside the
*.cu files -- we only mention some:

[1] dlpolycuda_twobody_max_device_iterations (default: 30000)
    Increasing this to match the system size may help.

[2] dlpolycuda_twobody_k2_unroll_slices (default: 6)
    Increasing the value, at the cost of gpu memory, may increase
    performance of the two-body force field functions
    considerably.  Available values 1-8, 10, 12 & 14.

[3] dlpolycuda_constraints_shake_ntcons_threshold (default: 0)
    The SHAKE algorithm will not be accelerated unless 'ntcons'
    is found >= than this value.

Please NOTE that the dlpolycuda_ options are used for development
purposes and changing their defaults may lead to undesired and
unforseen behaviours of the builds.

Contacts at ICHEC & STFC Daresbury Laboratory:
----------------------------------------------

Dr. R. Nestor     :: ruairi.nestor@ichec.ie
Dr. I.T. Todorov  :: ilian.todorov@stfc.ac.uk
