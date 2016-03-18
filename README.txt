DL_POLY_4.08
============

The source is in fully self-contained free formatted FORTRAN90+MPI2
code (specifically FORTRAN90 + TR15581 + MPI1 + MPI-I/O only).  Any
available functionality that dependent on external libraries; such as
NetCDF, PLUMED and OpenKIM; will require the user to satisfy any the
dependencies upon compilation and thus the use of `cmake' (refer to
README.md).  The non-extended, vanilla code complies with the NAGWare
and FORCHECK F90.

This version supports ALL features that are available in the
standard DL_POLY_Classic version with the exceptions of:
  (1) Rigid bodies (RBs) linked via (i) constraint bonds (CBs) or
      (ii) potential of mean field constraints (PMFs).
  (2) Truncated octahedral (imcon = 4), Rhombic Dodecahedral
      (imcon = 5) and Hexagonal Prism (imcon = 7) periodic boundary
      MD cell conventions (PBC).  Note that the last one is easily
      convertible to orthorhombic (imcon = 2), see utility/nfold.f90.
  (3) Classic Ewald and Hautman-Klein Ewald Coulomb evaluations.
  (4) Hyper-Dynamics: Temperature Accelerated Dynamics, Biased
      Potential Dynamics; Meta-Dynamics and solvation features.

No previous DL_POLY_3/4 feature is deprecated.  ALL NEW features are
documented in the "DL_POLY_4 User Manual".

References:
-----------
Thank you for using the DL_POLY_4 package in your work.  Please,
acknowledge our efforts by including the following references when
publishing data obtained using DL_POLY_4
  (1) "I.T. Todorov, W. Smith, K. Trachenko & M.T. Dove,
       J. Mater. Chem., 16 (20), 1911-1918 (2006)"
  (2) "I.J. Bush, I.T. Todorov & W. Smith,
       Comp. Phys. Commun., 175, 323-329 (2006)"
  (3) "H.A. Boateng & I.T. Todorov,
       J. Chem. Phys., 142, 034117 (2015)".

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
  (3) The `rpad' & `no strict' CONTROL options should be used with
      care especially in conjunction with the `variable timestep`
      option when iterative algorithms are present in the simulation.
      Such may be driven by a combination of options such as:
      `minimise', `ensemble npt', `ensemble nst' in the presence of
      constraints (associated with the `tolerance' and `mxshake'
      options) and/or core shells units (dealt by the relaxed shell
      model and associated with the `rlxtol' option) in the model
      system as defined in the FIELD file.
  (4) System stability is more easily compromised than the one when
      running DL_POLY_Classic!  When starting a new system it may be
      beneficial and considerate to make use of options that aim to
      provide extra safety; such as `l_scr', `scale', `cap forces'
      (`zero') and other tolerance related ones, and especially
      when in equilibration mode the `variable timestep` one in
      conjunction with an ensemble in Berendsen formulation!
  (5) The DL_POLY_4 parallel performance and efficiency are considered
      very-good-to-excellent as long as (i) all CPU cores are loaded
      with no less than 500 particles each and (ii) the major linked
      cells algorithm has no dimension less than 3.
  (6) Although DL_POLY_4 can be compiled in a serial mode, users are
      advised to consider DL_POLY_Classic as a suitable alternative to
      DL_POLY_4 when simulations are likely to be serial jobs for
      systems containing < 200 particles-per-processor.  In such
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
treat CB/PMF/RB entities in a more precise (symplectic) manner than
the LFV ones and thus not only have better numerical stability but
also produce more accurate dynamics.

Makefiles for legacy compilation:
--------------------------------
From within the `source' directory the user may compile the code by
selecting the appropriate Makefile from the `build' directory and
linking (or copying) it across by typing at the command line:

	ln -s ../build/Makefile_MPI Makefile

- intended for parallel execution on multi-processor platforms (an MPI
implementation is needed), or

	ln -s ../build/Makefile_SRL2 Makefile

- intended for serial execution (no MPI required)

followed by <Enter>.

If there is an `entry' in the Makefile for the particular combination
of architecture, compiler & MPI implementation, then the user may
instantiate the compilation by issuing at the command line:

	make `entry' BUILDER='Your Name Here'

and then pressing <Enter>.

Usually, the one named `hpc' is suitable for the majority of platforms.
Specifying your name as the person building the code is optional but
it may be useful metadata in the future for citation/acknowledgment
purposes when the results produced with the particular executable
are examined as part of a research programme.

In case the compilation process fails with a message about undefined
MPI related definitions in comms_module than it means that the MPI
libraries used in the process are old.  To rectify the matters you
will need to add `-DOLDMPI' to the FCFLAGS string within the relevant
Makefile `entry' before you try to compile again.

To find out the keywords for all available entries within the Makefile
issue:

	make

press <Enter> and then examine the Makefile entries corresponding to
the keywords reported.  If there is not a suitable entry then you
should seek advice from a computer scientist or the support staff of
the particular machine (HPC service).

The necessary components for the source compilation are:
  (1) a FORTRAN90 + TR15581 compliant compiler, i.e. a gfortran
      (requires a gcc version 4.2.0 or higher).  If the full PATH to
      it is not passed to the DEFAULT ENVIRONMENT PATH, then it MUST
      be explicitly supplied in the Makefile!
  (2) an MPI2 (or MPI1 + MPI-I/O) implementation, COMPILED for the
      architecture/OS and the targeted compiler.  Usually, this is
      encapsulated and badged as an mpif90 entry (or ftn on ARCHER)
      after appropriate module load on an HPC architecture.
      Otherwise, if the full PATH to these components is not passed
      to the DEFAULT ENVIRONMENT PATH then it MUST be explicitly
      supplied in the Makefile!
  (3) a MAKE command (Makefile interpreter in the system SHELL).
Note that (2) is not necessary for compilation in SERIAL mode!

By default, if the compilation process is successful then an executable
(build) will be placed in `../execute' directory (at the same level as
the `source' directory where the code is compiled).  Should the
`../execute' directory not exist then it will be created automatically
by the Makefile script.  The build may then be moved, renamed, etc. and
used as the user wishes.  However, when executed, the program will look
for input files in the directory of execution!

Compilation on Windows:
-----------------------
The best way to get around it is to install cygwin on the system
(http://www.cygwin.com/) to emulate a UNIX/Linux like environment
and then use the "make" command.  During the cygwin installation please
make sure that the "make" and "gfortran" components are specifically
opted for components (as they may not be included as default ones) in
the install.  A potential problem for Windows based FORTRAN compilers,
one may encounter, is that the compiler may not pick symbolic links.
This can be resolved by substituting the soft links with hard in the
Makefile.  For parallel compilation all "openMPI" components must also
be opted for during the cygwin install!!!

Contacts at STFC Daresbury Laboratory
-------------------------------------
Dr. I.T. Todorov :: ilian.todorov@stfc.ac.uk
