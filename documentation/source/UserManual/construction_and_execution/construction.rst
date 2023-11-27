Constructing : an Overview
==========================

Constructing the Standard Versions
----------------------------------

DL_POLY_4 was designed as a package of useful subroutines rather than a
single program, which means that users are to be able to construct a
working simulation program of their own design from the subroutines
available, which is capable of performing a specific simulation.
However, we recognise that many, perhaps most, users will be content
with creating a standard version that covers all of the possible
applications with DL_POLY_4 native functionalities and may be a few
external ones for this reason we have only provided the necessary tools
to assemble such a version. The methods of creating the standard
versions is described in detail in this chapter, however a brief
step-by-step description follows.

#. DL_POLY_4 is supplied as a ZIP compressed file. This must
   uncompressed to create the DL_POLY_4 directory
   (Section `[directory-structure] <#directory-structure>`__).

#. The next step is to compile DL_POLY_4 either by using the traditional
   Makefiles or by creating OS customised ones using **cmake**. In
   either case **make** will be used to produce a binary executable (see
   Section :ref:`Compiling and Running<compilation>`), which as a default will
   be named DLPOLY.Z and located in the *execute* subdirectory if
   compiling in traditional mode or in the place of using **cmake**.

#. DL_POLY_4 also has a Java GUI. The files for this are stored in the
   subdirectory *java*. Compilation of this is simple and requires
   running the javac compiler and the jar utility. Details for these
   procedures are provided in the GUI manual
   :cite:`smith-gui`.

#. To run the executable for the first time you require the files
   CONTROL, FIELD and CONFIG; and possibly a few tabulated files (TABLE,
   TABEAM, etc.). These must be present in the directory from which the
   program is executed. (See Section :ref:`The INPUT Files<input-files>`
   for the description of the input files.)

#. Executing the program will produce the files OUTPUT, STATIS, REVCON
   and REVIVE; and optionally a constellation of others in the executing
   directory depending on user options in CONTROL (HISTORY, RDFDAT,
   ZDNDAT, MSDTMP, REFERENCE, DEFECTS, etc.). (See
   Section :ref:`The OUTPUT Files<output-files>` for the description of
   the output files.)

This simple procedure is enough to create a standard version to run most
simulations. There may however be some difficulty with array sizes.
DL_POLY_4 contains features which allocate arrays after scanning the
input files for a simulation. Sometimes these initial estimates are
insufficient for a long simulation when, for example, the system volume
changes markedly during the simulation or when a system is artificially
constructed to have a non-uniform density. Usually, simply restarting
the program will cure the problem, but sometimes, especially when the
local atom density is somewhat higher than the global one or the system
undergoes some form of clustering and the distribution of bonded-like
interactions is far from uniform, it may be necessary to amend the array
sizes in accordance with the error message obtained. To trigger
lengthening of the density dependent global arrays the user may use the
**densvar** option in the CONTROL
(Section :ref:`The CONTROL File<control-file>`) file. However, lengthening
these arrays will require a larger amount of memory from the execution
machine for the simulation, which it may not be able to provide. See
Section :ref:`Source Code<file-structure>` for more insight on the
DL_POLY_4 source code structure.

Constructing Non-standard Versions
----------------------------------

In constructing a non-standard DL_POLY_4 simulation program, the first
requirement is for the user to write a program to function as the root
segment. The root segment ``/VV/dl_poly`` is placed in the *source*
directory and contains the set-up and close-down calls for a molecular
dynamics simulation. It is the routine that first opens the OUTPUT file
(Section :ref:`The OUTPUT Files<output-files>`), which provides the summary
of the job. The root program calls the “molecular dynamics cycle”
routines implementing the VV. These routines contain major routines
required to perform the simulation, control the normal “molecular
dynamics cycle” and monitor the *cpu* and *memory* usage. They also
bring about a controlled termination of the program if the *cpu* usage
approaches the allotted job time within a pre-set closure time and/or if
the *memory* usage approaches the allocated limit for density dependent
arrays. Users are recommended to study the aforementioned root
directories as a model for other implementations of the package they may
wish to construct. Some advise on hierarchies of all the DL_POLY_4
subroutines can be found in
Section :ref:`File Structure<file-structure>`.

Should additional functionality be added to DL_POLY_4 by the user, the
``set_bounds`` routine (and its support subroutines) may need modifying
to allow specification of the dimensions of any new arrays.

Any molecular dynamics simulation performs five different kinds of
operation: initialisation; forces calculation; integration of the
equations of motion; calculation of system properties; and job
termination. It is worth considering these operations in turn and to
indicate which DL_POLY_4 routines are available to perform them. We do
not give a detailed description, but provide only a guide. Readers are
recommended to examine the different routines described in the DL_POLY_4
User Manual for further details (particularly regarding further
dependencies i.e. additional routines that may be called).

The following outline assumes a system containing flexible molecules
held together by rigid bonds.

Initialisation requires firstly that the program determine what platform
resources are made available to the specific simulation job. This is
done by the DL_POLY_4 routine ``map_domains`` in ``domains_module`` that
attempts to allocate and map the resources (nodes in parallel) in
compliance with the DD :index:`strategy<parallelisation>`. ``map_domains`` is called within the
routine ``set_bounds``, which also sets the necessary limits for various
simulation array sizes and all global variables as declared in
``setup_module`` to convenient values based on a rough scan through the
CONFIG, CONTROL, FIELD and optionally TABLE and TABEAM
(Section :ref:`The INPUT Files<input-files>`) files. The routine also calls
the ``read_config`` routine to obtain atomic positions and optionally
velocities and forces from the CONFIG file. After allocation of all
necessary simulation arrays and variables (with compulsory
initialisation to “zero” value), the job control information is
required; this is obtained by the routine ``read_control``, which reads
the CONTROL file. The description of the system to be simulated – the
types of atoms and molecules present and the intermolecular forces – are
obtained by the ``read_field`` routine, which reads the FIELD file. The
``system_init`` routine is called next to initialise various simulation
arrays and variables with the data available so far and detects if the
job is a restart of previous simulation run. If so it reads the REVOLD
(Section :ref:`The REVOLD File<revold-file>`) to supply some arrays and
variables with the necessary values as saved from the previous job. The
domain halo is constructed immediately afterwards by the routine
``set_halo_particles``. After gathering all these data, bookkeeping and
exclusion arrays are created for the intramolecular and site related
interactions (core-shell, constraint and tether units) by the
``build_book_intra`` and ``build_excl_intra`` routines. Lastly, the
thermodynamic properties of the system are checked and set by the
``set_temperature`` routine (which also generates the initial velocities
if required to do so).

The calculation of the pair-like forces is carried out in the
``two_body_forces`` routine and represents the main part of any
simulation. For calculation of the two-body contributions to the atomic
forces, the :index:`Verlet neighbour list` is constructed by the
``link_cell_pairs`` routine using link-cell lists. Special measures are
taken so that the list excludes: (i) pairs of atoms that are both in a
*frozen state* as well as (ii) pairs in which one of the atoms has the
other in its *exclusion list*. The last is built by ``build_excl_intra``
where the specifications of bond-like interactions in the FIELD file are
processed. Various other subroutines are then called to calculate
specific contributions by different interactions. For example;
``vdw_forces`` for the short-range (van der Waals) 
:index:`forces<potential;van der Waals>`
(Section :ref:`Short Ranged (van der Waals) Potentials<vdw>`), ``metal_lrc``, ``metal_ld_compute`` and
``metal_forces`` for the metal :index:`interactions<potential;metal>`
(Section :ref:`Metal Potentials<metal>`), and ``ewald_spme_forces``,
``ewald_real_forces``, ``ewald_frzn_forces`` and ``ewald_excl_forces``
for the Coulombic forces (Section :ref:`Long Ranged Electrostatic (coulombic) Potentials<coulomb>`).

Higher order :index:`intermolecular<potential;intermolecular>`, 
site-related and :index:`intramolecular<potential;intramolecular>` forces
require the routines
``tersoff_forces``, ``three_body_forces``, ``four_body_forces``,
``core_shell_forces`` or ``core_shell_relax``, ``tethers_forces``,
``bonds_forces``, ``angles_forces``, ``dihedrals_forces`` and ``inversions_forces``.
The routines
``external_field_apply`` and ``external_field_correct`` are required
if the simulated system has an external force field (e.g.
electrostatic field) operating.

To help with equilibration simulations, routines such as ``cap_forces``,
``zero_k_optimise`` and ``minimise_relax`` are sometimes required to
reduce the magnitude of badly equilibrated forces and to steer the MD
system towards an equilibrium state.

Integration of the equations of motion is handled by one of the routines
listed and described in
Chapter :ref:`Integration Algorithms<integration-algorithms>`.

As mentioned elsewhere, DL_POLY_4 does not contain many routines for
computing system properties during a simulation. Radial distributions
may be calculated however, by using the routines ``rdf_collect``,
``rdf_excl_collect``, ``rdf_frzn_collect`` a and ``rdf_compute``.
Similarly, Z-density distributions may be calculated by using the
routines ``z_density_collect`` and ``z_density_compute``, while velocity
autocorrelation functions may be calculated using the routines
``vaf_collect`` and ``vaf_compute``. Ordinary thermodynamic quantities
are calculated by the routine ``statistics_collect``, which also writes
the STATIS file (Section :ref:`The STATIS File<statis-file>`). Routine
``trajectory_write`` writes the HISTORY
(Section :ref:`The HISTORY<history-file>`) file for later (postmortem)
analysis. Routine ``defects_write`` writes the DEFECTS
(Section :ref:`The DEFECTS File<defects-file>`) file for later (postmortem)
analysis. Routine ``msd_write`` writes the MSDTMP
(Section :ref:`The MSDTMP File<msdtmp-file>`) file for later (postmortem)
analysis. Routine ``rsd_write`` writes the RSDDAT
(Section :ref:`The RSDDAT File<rsddat-file>`) file for later (postmortem)
analysis.

Job termination is handled by the routine ``statistics_result`` which
writes the final summaries in the OUTPUT file and dumps the restart
files REVIVE and REVCON (Sections :ref:`The REVIVE File<revive-file>` and
`The REVCON File<revcon-file>` respectively).
