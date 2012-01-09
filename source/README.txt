DL_POLY_4.03
============

The source is in fully self-contained free formatted FORTRAN90+MPI2
code (specifically FORTRAN90 + TR15581 + MPI1 + MPI-I/O only).  The
available NetCDF functionality makes the extended code dependent upon
it.  The non-extended code complies with the NAGWare f95 and FORCHECK
f90 standards with exception of the FORTRAN2003 feature TR15581, which
is very rarely unavailable in the nowadays FORTRAN95 compilers.

This version supports ALL features that are available in the
standard DL_POLY_Classic version with the exceptions of:
  (1) RIDGID BODIES linked by constraint bonds (CB) or
      potential of mean field (PMF) constraints.
  (2) Truncated octahedral (imcon = 4), Rhombic Dodecahedral
      (imcon = 5) and Hexagonal Prism (imcon = 7) periodic boundary
      conventions.
  (3) Classic Ewald and Hautman-Klein Ewald Coulomb evaluations.
  (4) Temperature Accelerated Dynamics, Hyper-Dynamics and
      solvation energies.

No previous DL_POLY_3/4 feature is deprecated.  ALL NEW features are
documented in the "DL_POLY_4 User Manual".

Reference:
---------
Thank you for using the DL_POLY_4 package in your work.  Please,
acknowledge our efforts by including the following reference when
publishing data obtained using DL_POLY_4: "I.T. Todorov, W. Smith,
K. Trachenko & M.T. Dove, J. Mater. Chem., 16, 1611-1618 (2006)".

Warnings:
---------
  (1) DL_POLY_4 can produce index ordered REVCON, HISTORY and MSDTMP
      files which are restartable by DL_POLY_Classic.  Although such
      printed outputs look unscrambled, the actual printing process
      is not.  Unscrambled printing is slightly more expensive than
      natural (scrambled) printing.  The cost time-wise is little,
      < 1%, but HD space-wise is approximately 20%.  This is due to
      the necessary addition of blanks at the end of data record,
      included to align the (ASCII) lines of output files (human
      readable) to a constant length.  Printing scrambled outputs
      is optional.  Note that these too have blanks aligned records.
      The parallel I/O ensures (i) writing speeds of 10^5 to 10^6
      particle per second with optimal number of writers and (ii)
      reading speeds of 10^4 to 10^5 particles per second per reader.
      For more information on I/O options consult the user manual.
  (2) REVIVE files produced by different versions are not compatible.
      Furthermore, restarting runs across different sub-versions
      may not be possible.
  (3) The DL_POLY_4 parallel performance and efficiency are considered
      very-good-to-excellent as long as (i) all CPU cores are loaded
      with no less than 500 particles each and (ii) the major linked
      cells algorithm has no dimension less than 4.
  (4) Although DL_POLY_4 can be compiled in a serial mode, users are
      advised to consider DL_POLY_Classic as a suitable alternative to
      DL_POLY_4 when simulations are likely to be serial jobs for
      systems containing < 500 particles-per-processor.  In such
      circumstances, with both codes compiled in serial mode, the
      difference in performance, measured by the time-per-timestep
      ratio [DL_POLY_Classic(t)-DL_POLY_4(t)]/DL_POLY_Classic(t),
      varies in the range -5:+5%.  This variation depends strongly on
      the system force-field complexity and very weakly on the system
      size.

Integration Defaults:
---------------------
The default ensemble is NVE.

The default integration scheme is Trotter derived Velocity Verlet
(VV), although Leapfrog Verlet (LFV) is also available.  VV is
considered superior (to LFV) since:
  (1) Integration can be developed in symplectic manner for certain
      ensembles, such as: NVE, NVEk (NVT Evans) as well as all
      Nose-Hoover ensembles (NVT, & NPT & NsT when there is no
      external field applied on the system, otherwise they do not
      conserve the phase space volume) and MTK ensembles (NPT & NsT).
  (2) All ensemble variables are updated synchronously and
      thermodynamic quantities and estimators are exact at the every
      step, whereas in LFV particle velocities and thermostat and
      barostat friction velocities are half an integration time-step
      behind the rest of the ensemble variables and due to this
      certain estimators are approximated at full timestep.
  (3) It offers better numerical stability and faster convergence
      when (i) constraint solvers (CB/PMF: RATTLE/VV versus SHAKE/LFV)
      are involved and/or (ii) RB dynamics is integrated.

The LFV integration may take less cpu time than the VV one for the
certain ensembles - type of system (CB/PMF/RB) and type of ensemble
dependent.  Usually, LFV is slightly faster than VV when CB/PMF/RB
are present in the system.  The relative performance between the LVF
and VV integration (per timestep) is observed to vary in the limits
*** [LFV(t)-VV(t)]/VV(t) = -5:+5% ***.  However, the VV algorithms
treat CB/PMF/RB entities in more precise (symplectic) manner than
the LFV ones and thus not only have better numerical stability but
also produce more accurate dynamics.

Makefiles:
----------
From within the `source' directory the user may compile the code by
selecting the appropriate Makefile from the `build' directory:

"cp ../build/Makefile_MPI Makefile"  (for parallel execution MPI is needed)

or

"cp ../build/Makefile_SRLx Makefile" (for serial execution - no MPI needed)

Note that in `comms_module.f90' it is crucial that line 13 reads as:
`Use mpi_module'     for serial compilation and
`Use mpi'            for parallel compilation (which is the default)

If the parallel OS environment, you are compiling on, is not fully F90
compatible then the `Use mpi' entry in `comms_module.f90' will be
interpreted as erroneous.  This is easily overcome by commenting out
`Use mpi' and inserting "Include 'mpif.h'" after `Implicit None'.

If there is an `entry' in the Makefile for the particular combination
of architecture, compiler & MPI library, then the user may
instantiate the compilation by:

"make `entry'"

If there is not a suitable entry, the user should advise with a
computer scientist or the administrator of the particular machine.

The necessary components for the source compilation are:
  (1) a FORTRAN90 compliant compiler (if the full PATH to it is not
      passed to the DEFAULT ENVIRONMENT PATH, then it MUST be
      explicitly supplied in the Makefile)
  (2) MPI2 (or MPI1 + MPI-I/O) libraries COMPILED for the architecture
      and the targeted compiler (if the full PATH to these is not
      passed to the DEFAULT ENVIRONMENT PATH, then it MUST be
      explicitly supplied in the Makefile)
  (3) a MAKE command (Makefile interpreter in the system SHELL)

Note that (2) is not necessary for compilation in SERIAL mode!

By default, if compilation is successful, an executable (build) will
be placed in "../execute" directory (at the same level as the
directory where the code is compiled).  Should it not exist one will
be created automatically.  The build can then be moved, renamed, etc.
and used as the user wishes.  However, when executed, the program will
look for input files in the directory of execution!

Serial Compilation on Windows:
------------------------------
The best way to get around it is to install cygwin on the system
(http://www.cygwin.com/) to emulate a UNIX/Linux like environment and
then use the "make" command.  During cygwin installation make sure that
make and gfortran are included in the install.  A potential problem
for Windows based FORTRAN compilers, you may encounter, is that the
compiler may not pick symbolic links.  To resolve this, you will
have to use hard linking in the Makefile.

Compiling with NetCDF functionality:
------------------------------------
The targeted Makefile needs the following substitution within before
attempting compilation:

"netcdf_modul~.o -> netcdf_module.o"

Note that suitable entry may need to be created within the Makefile
so that it matches the particular combination of architecture,
compiler, MPI library & netCDF library.

Contacts at STFC Daresbury Laboratory:
--------------------------------------
Dr. I.T. Todorov :: ilian.todorov@stfc.ac.uk
