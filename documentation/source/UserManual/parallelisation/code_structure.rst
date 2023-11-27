.. _source-code_sec:

Source Code
===========

.. _modularisation:

Modularisation Principles
-------------------------

Modules in DL_POLY_4 are constructed to define parameters and variables
(scalars and arrays) and/or develop methods that share much in common.
The division is far from arbitrary and module interdependence is reduced
to minimum. However, some dependencies exist which leads to the
following division by groups in hierarchical order:

-  **precision module**:: ``kinds_f90``

   The precision module defines the working precision ``wp`` of all real
   variables and parameters in . By default it is set to 64-bit (double)
   precision. If the precision is changed, the user must check whether
   the specific platform supports it and make sure it is allowed for in
   the MPI implementation. If all is OK then the code must be
   recompiled.

-  **MPI module**:: ``mpi_module``

   The MPI module implements all MPI functional calls used in . It is
   only used when DL_POLY_4 is to be compiled in serial mode.

-  **communication module**:: ``comms_module`` (``mpi_module``)

   The communication module defines MPI related parameters and develops
   MPI related functions and subroutines such as: initialisation and
   exit; global synchronisation, sum, maximum and minimum; node ID and
   number of nodes; simulation time. It is dependent on ``kinds_f90``
   and on ``mpi_module`` if MPI is emulated for DL_POLY_4 compilation in
   serial mode. The ``mpi_module`` implements all MPI functional calls
   used in .

-  **global parameters module**:: ``setup_module``

   The global parameters module holds all important global variables and
   parameters (see above). It is dependent on ``kinds_f90``.

-  **parse module**:: ``parse_module``

   The parse module develops several methods used to deal with textual
   input: ``get_line strip_blanks lower_case get_word word_2_real``.
   Depending on the method dependencies on
   ``kinds_f90 comms_module setup_module domains_module`` are found.

-  **development module**:: ``development_module``

   The development module contains several methods used to help with
   testing and debugging . Depending on the method dependencies on
   ``kinds_f90 comms_module setup_module domains_module`` are found.

-  **netCDF module**:: ``netcdf_module``

   The netCDF module contains all important netCDF functions and global
   variables in DL_POLY_4 context. It is dependent on ``kinds_f90``.

-  **I/O module**:: ``io_module``

   The I/O module contains all important global variables that define
   the I/O methods and types used in the package and contains basic
   routines essential for the I/O in . It is dependent on ``kinds_f90``.

-  **domains module**:: ``domains_module``

   The domains module defines DD parameters and maps the available
   computer resources on a DD grid. The module does not depend on
   previous modules but its mapping subroutine is dependent on
   ``kinds_f90`` and ``comms_module``.

-  **site module**:: ``site_module``

   The site module defines all site related arrays (FIELD) and is
   dependent on ``kinds_f90`` only. However, it also develops an
   allocation method that is dependent on ``setup_module``.

-  **configuration module**:: ``config_module``

   The configuration module defines all configuration related arrays
   (CONFIG) and is dependent on ``kinds_f90`` only. However, it also
   develops an allocation method that is dependent on ``setup_module``.

-  **vnl module**:: ``vnl_module``

   The Verlet neighbour list (VNL) module defines all VNL related
   control variables and arrays needed for the VNL conditional update
   functionality, and is dependent on ``kinds_f90`` only. However, it is
   assisted by a ``vnl_check`` routine that is dependent on more
   modules.

-  **defects modules**:: ``defects_module defects1_module``

   The defects modules define all defects and configuration related
   arrays (REFERENCE) and are dependent on ``kinds_f90`` only. However,
   they also develop allocation methods that are dependent on
   ``setup_module``.

-  **dpd module**:: ``dpd_module``

   The Dissipative Particle Dynamics (DPD) module defines all DPD
   related control variables and arrays needed for the DPD
   functionality, and is dependent on ``kinds_f90`` only. However, it
   also develops an allocation method that is dependent on
   ``setup_module``.

-  **electrostatic modules**:: ``ewald_module mpoles_module``

   This modules define all variables and arrays needed for the SPME (i)
   refreshment k-space driven properties (ii) and multipola
   relectrostatics control variable and arrays in the DL_POLY_4 scope
   when. They depend on ``kinds_f90`` and but their allocation methods
   on ``setup_module``.

-  **inter-molecular interactions modules**:: vdw_module metal_module
   tersoff_module three_body_module four_body_module

   The intermolecular modules define all variables and potential arrays
   needed for the calculation of the particular interaction in the
   DL_POLY_4 scope. They depend on ``kinds_f90``. Their allocation
   methods depend on ``setup_module``.

-  **extra-molecular interactions modules**:: kim_module plumed_module

   These modules define all variables, arrays and functions needed
   OpenKIM and PLUMED integrable into DL_POLY_4 plugins. They depend on
   ``kinds_f90``. Their allocation methods depend on ``setup_module``.

-  **intra-molecular interactions and site-related modules**::
   ``rdf_module z_density_module core_shell_module constraints_module pmf_module rigid_bodies_module tethers_module bonds_module angles_module dihedrals_module inversions_module``

   These modules define all variables and potential or statistical grid
   arrays needed for the calculation of the particular interaction or
   distribution function in the DL_POLY_4 scope. They all depend on
   ``kinds_f90`` with allocation methods depending on ``setup_module``.

-  **external field module**:: ``external_field_module``

   This module defines all variables and potential arrays needed for the
   application of an external field in the DL_POLY_4 scope. It depends
   on ``kinds_f90`` and its allocation method on ``setup_module``.

-  **langevin module**:: ``langevin_module``

   This module defines all variables and arrays needed for the
   application of NPT and N\ :math:`\underline{\underline{\mathbf{\sigma}}}`\ T Langevin routines
   in the DL_POLY_4 scope. It depends on ``kinds_f90`` and its
   allocation method on ``setup_module``.

-  **minimise module**:: ``minimise_module``

   This module defines all variables and arrays needed for the
   application of a Conjugate Gradient Method minimisation routine in
   the DL_POLY_4 scope. It depends on ``kinds_f90`` and its allocation
   method on ``setup_module``.

-  **msd module**:: ``msd_module``

   This module globalises a CONTROL variable.

-  **statistics module**:: ``statistics_module``

   This module defines all variables and arrays needed for the
   statistical accountancy of a simulation in . It depends on
   ``kinds_f90`` and its allocation methods on ``setup_module`` and
   ``comms_module``.

-  **greenkubo module**:: ``greenkubo_module``

   This module defines all variables and arrays needed for calculation
   of Green-Kubo relations during a simulation in . It depends on
   ``kinds_f90`` and its allocation methods on ``setup_module``.

-  **kinetic module**:: ``kinetic_module``

   The kinetic module contains a collection of routines for the
   calculation of various kinetic properties. It is dependent on
   ``kinds_f90``.

-  **DaFT module**:: ``gpfa_module parallel_fft``

   These modules contain all necessary functionality for DL_POLY_4 DaFT
   and it GPFA 1D FFT dependence. They have dependencies on
   ``kinds_f90``, ``comms_module.f90`` and ``setup_module.f90``.

.. _file-structure:

File Structure
--------------

Generally, the DL_POLY_4 file structure can be divided into four groups
as follows:

-  **general files** in the *source* directory

-  **SERIAL specific** files in the *source/SERIAL* directory

The files in each group are listed in hierarchal order as closely as
possible as examplified in the relevant DL_POLY_4 Makefies in the
*build* subdirectory. The further down the category the file, the more
dependent it is on the files listed above it.

Module Files
------------

The DL_POLY_4 module files contain all global variables (scalars and
arrays) and parameters as well as some general methods and generic
functions intrinsically related to the purpose or/and contents of the
specific module. The file-names and the methods or/and functions
developed in them have self-explanatory names. More information of their
purpose can be found in their headers.

The rest of files in DL_POLY_4 are dependent on the module files in
various ways. The dependency relation to a module file is explicitly
stated in the declaration part of the code.

General Files
-------------

The DL_POLY_4 general files are common to both MPI and SERIAL version of
the code. In most cases, they have self-explanatory names as their order
is matched as closely as possible to that occurring in the main segment
of the code - ``dl_poly``. Only the first five files are exception of
that rule; ``warning`` and ``error`` are important reporting subroutines
that have call points at various places in the code, and
``numeric_container``, and ``spme_container`` are containers of simple
functions and subroutines related in some way to their purpose in the
code.

SERIAL Specific Files
---------------------

These implement an emulation of some general MPI calls used in source
code when compiling in serial mode as well as some modified counterparts
of the general files changed to allow for faster and/or better memory
optimised serial execution. Names are self-explanatory.

Comments on MPI Handling
------------------------

Only a few files make explicit calls to MPI routines.

.. _parameters:

Comments on ``setup_module``
----------------------------

The most important module, by far, is ``setup_module``, which holds the
most important global parameters and variables (some of which serve as
“parameters” for global array bounds, set in set_bounds). A brief
account of these is given below:

.. list-table::
  :header-rows: 1 

  * - Parameter 
    - Value 
    - Function 
  * - ``DLP_VERSION``
    - string
    -  version string - number 
  * - ``DLP_RELEASE`` 
    - string 
    - release string - date 
  * - ``pi`` 
    - 3.14159265358979312 
    - :math:`\pi` constant
  * - ``twopi`` 
    - 6.28318530717958623 
    - :math:`2 \pi` constant
  * - ``fourpi`` 
    - 12.56637061435917246 
    - :math:`4 \pi` constant 
  * - ``sqrpi`` 
    - 1.772453850905588 
    - :math:`\sqrt[2]{\pi}` constant
  * - ``rtwopi`` 
    - 0.15915494309189535
    - :math:`\frac{1}{2 \pi}` constant
  * - ``rt2``
    - 1.41421356237309515
    - :math:`\sqrt[2]{2}` constant
  * - ``rt3`` 
    - 1.73205080756887719 
    - :math:`\sqrt[2]{3}` constant 
  * - ``r4pie0``
    - 138935.4835
    - electrostatics conversion factor to internal units, i.e. :math:`\frac{1}{4 \pi \epsilon_{o}}`
  * - ``boltz``
    - 0.831451115
    - Boltzmann constant in internal units also used as Kelvin/Boltzmann energy unit (very rarely used)
  * - ``engunit`` 
    - *variable* 
    - the system energy unit
  * - ``eu_ev`` 
    - 9648.530821
    - eV energy unit (most used)
  * - ``eu_kcpm``
    - 418.4 kcal/mol
    - energy unit (often used)
  * - ``eu_kjpm``
    - 100.0 kJouls/mol
    - energy unit (rarely used)
  * - ``prsunt`` 
    - 0.163882576 
    - conversion factor for pressure from internal units to katms
  * - ``tenunt`` 
    - 1.660540200 
    - conversion factor for surface tension from internal units to dyn/cm
  * - ``del_max``
    - 0.01 
    - maximum bin sizes in Angstroms for distance grids
  * - ``delth_max``
    - 0.20
    - maximum bin sizes in degrees for angle grids
  * - ``nread`` 
    - 5
    - main input channel
  * - ``nconf``  
    - 11 
    - configuration file input channel
  * - ``nfield``
    - 12 
    - force field input channel
  * - ``ntable`` 
    - 13 
    - tabulated potentials file input channel
  * - ``nrefdt`` 
    - 14
    - reference configuration input channel
  * - ``nrite`` 
    - 6 
    - main output channel 
  * - ``nstats`` 
    - 21 
    - statistical data file output channel 
  * - ``nrest`` 
    - 22 
    - output channel accumulators restart dump file 
  * - ``nhist`` 
    - 23 
    - trajectory history file channel
  * - ``ndefdt`` 
    - 24 
    - output channel for defects data file 
  * - ``nrdfdt`` 
    - 25 
    - output channel for RDF data
  * - ``nzdfdt`` 
    - 26 
    - output channel for Z-density data file
  * - ``nrsddt`` 
    - 27 
    - output channel for displacements data files
  * - ``npdfdt`` 
    - 28 
    - output channel for raw PDF files 
  * - ``ngdfdt`` 
    - 29 
    - output channel for normalised RDF data files
  * - ``nvafdt`` 
    - 28
    -  output channel for VAF files
  * - ``nmpldt``
    - 29 
    - output channel for the PLOLES data file
  * - ``seed(1:3)`` 
    - *variable* 
    - pair of seeds for the random number generator 
  * - ``lseed`` 
    - *variable* 
    - logical swich on/off indicator for seeding 
  * - ``mxsite`` 
    - *variable*
    - max number of molecular sites
  * -  ``mxatyp``
    - *variable*
    - max number of unique atomic types
  * - ``mxtmls`` 
    - *variable* 
    - max number of unique molecule types
  * - ``mxexcl`` 
    - *variable*
    - max number of excluded interactions per atom
  * - ``mxompl``
    - *variable*
    - max number of multipolar order specified
  * - ``mximpl`` 
    - *variable* 
    - max number of multipolar total momenta for this order
  * - ``mxspl`` 
    - *variable*  
    - SPME FFT B-spline order
  * - ``mxspl1`` 
    - *variable* 
    - SPME FFT B-spline possible extension when :math:`r_{\rm pad}>0`
  * - ``kmaxa`` 
    - *variable* 
    - SPME FFT amended array dimension (a direction)
  * - ``kmaxb`` 
    - *variable* 
    - SPME FFT amended array dimension (b direction)
  * - ``kmaxc`` 
    - *variable*
    - SPME FFT amended array dimension (c direction)
  * - ``kmaxa1`` 
    - *variable*
    - SPME FFT original array dimension (a direction)
  * - ``kmaxb1``
    - *variable* 
    - SPME FFT original array dimension (b direction)
  * - ``kmaxc1``
    - *variable*
    - SPME FFT original array dimension (c direction)
  * - ``mxtshl``
    - *variable* 
    - max number of specified core-shell unit types in system
  * - ``mxshl`` 
    - *variable*
    - max number of core-shell units per node
  * - ``mxfshl``
    - *variable* 
    - max number of related core-shell units (1+1)
  * - ``mxtcon`` 
    - *variable*
    - max number of specified bond constraints in system
  * - ``mxcons`` 
    - *variable*
    - max number of constraint bonds per a node
  * - ``mxfcon`` 
    - *variable* 
    - max number of related constraint units (6+1)
  * - ``mxlshp``
    - *variable* 
    - max number of shared particles per node :math:`\texttt{Max} (2~\frac{\texttt{mxshl}}{2},2~\frac{\texttt{mxcons}}{2},\frac{\texttt{mxlrgd}~*~\texttt{mxrgd}}{2})`
  * - ``mxproc`` 
    - *variable* 
    - number of neighbour nodes in DD hypercube (26)
  * - ``mxtpmf(1:2)``
    - *variable*
    - max number of specified particles in a PMF unit (1:2)
  * - ``mxpmf``
    - *variable*
    - max number of PMF constraints per a node
  * - ``mxfpmf`` 
    - *variable*
    - max number of related PMF units (1+1)
  * - ``mxtrgd`` 
    - *variable*
    - max number of types RB units
  * - ``mxrgd`` 
    - *variable*
    - max number of RB units per node
  * - ``mxlrgd`` 
    - *variable*
    - max number of constituent particles of an RB unit
  * - ``mxfrgd`` 
    - *variable* 
    - max number of related RB units (1+1)
  * - ``mxtteth`` 
    - *variable*
    - max number of specified tethered potentials in system
  * - ``mxteth`` 
    - *variable* 
    - max number of tethered atoms per node
  * - ``mxftet`` 
    - *variable*
    -  max number of related tether units (1+1)
  * - ``mxpteth`` 
    - *variable*
    - max number of parameters for tethered potentials (3)
  * - ``mxtbnd`` 
    - *variable*
    - max number of specified chemical bond potentials in system
  * - ``mxbond``
    - *variable*
    - max number of chemical bonds per node
  * - ``mxfbnd``
    - *variable*
    - max number of related chemical bonds (1+(6*(6+1))/2)
  * - ``mxpbnd``
    - *variable* 
    - max number of parameters for chemical bond potentials (4)
  * - ``mxgbnd`` 
    - *variable*
    - max number of grid points in chemical bond pot. arrays (:math:`> 1004`)
  * - ``mxtang`` 
    - *variable*
    - max number of specified bond angle potentials in system
  * - ``mxangl``
    - *variable*
    - max number of bond angles per node
  * - ``mxfang`` 
    - *variable*
    - max number of related bond angles (1+(6*(6+1))/2)
  * - ``mxpang`` 
    - *variable*
    - max number of parameters for bond angle potentials (6)
  * - ``mxgang``
    - *variable*
    - max number of grid points in bond angle pot. arrays (:math:`> 1004`)
  * - ``mxtdih`` 
    - *variable*
    - max number of specified dihedral angle potentials in system
  * - ``mxdihd``
    - *variable*
    - max number of dihedral angles per node
  * - ``mxfdih``
    - *variable*
    - max number of related dihedral angles (1+((6-2)6*(6+1))/2)
  * - ``mxpdih``
    - *variable*
    - max number of parameters for dihedral angle potentials (7)
  * - ``mxgdih``  
    - *variable* 
    - max number of grid points in dihedral angle pot. arrays (:math:`> 1004`)
  * - ``mxtinv``
    - *variable*
    - max number of specified inversion angle potentials in system
  * - ``mxinv`` 
    - *variable* 
    - max number of inversion angles per node
  * - ``mxfinv`` 
    - *variable* 
    - max number of related inversion angless (1+(6*(6+1))/4)
  * - ``mxpinv`` 
    - *variable* 
    - max number of parameters for inversion angle potentials (3)
  * - ``mxginv``
    - *variable* 
    -  max number of grid points in inversion angle pot. arrays (:math:`> 1004`)
  * - ``mxrdf`` 
    - *variable*
    - max number of pairwise RDF in system
  * - ``mxgrdf`` 
    - *variable*
    - number of grid points for RDF and Z-density arrays (:math:`> 1004`)
  * - ``mxgele`` 
    - *variable* 
    - max number of grid points for ewald exclusion potential arrays
  * - ``mxgusr`` 
    - *variable*
    - number of grid points for umbrella sampling restraint’s RDF (:math:`> 1004`)
  * - ``mxvdw`` 
    - *variable*
    - max number of van der Waals potentials in system
  * - ``mxpvdw`` 
    - *variable* 
    - max number of van der Waals potential parameters (5)
  * - ``mxgvdw`` 
    - *variable* 
    -  max number of grid points in vdw potential arrays (:math:`> 1004`)
  * - ``mxmet`` 
    - *variable*
    - max number of metal potentials in system
  * - ``mxmed`` 
    - *variable*
    - max number of metal density potentials in system
  * - ``mxmds`` 
    - *variable* 
    - max number of metal extra density potentials in system
  * - ``mxpmet`` 
    - *variable*
    - max number of metal potential parameters (9)
  * - ``mxgmet``
    - *variable* 
    - max number of grid points in metal potential arrays (:math:`> 1004`)
  * - ``mxter`` 
    - *variable*
    - max number of Tersoff potentials in system
  * - ``mxpter`` 
    - *variable*
    -  max number of Tersoff potential parameters (11)
  * - ``mxgter`` 
    - *variable*
    - max number of grid points in tersoff potential arrays (:math:`> 1004`)
  * - ``mxgrid``
    - *variable* 
    - max number of grid points in potential arrays (:math:`> 1004`)
  * - ``mxtana``
    -  *variable* 
    - max number of PDFs per type
  * - ``mxgana`` 
    - *variable* 
    - max number of grid points for PDFs arrays
  * - ``mxgbnd``
    - *variable* 
    - max number of grid points for chemical bonds PDFs
  * - ``mxgang`` 
    - *variable*
    - max number of grid points for bond angles PDFs
  * - ``mxgdih`` 
    - *variable*
    - max number of grid points for dihedral angles PDFs
  * - ``mxginv`` 
    - *variable*
    - max number of grid points for inversion angles PDFs
  * - ``mxtbp``
    - *variable*
    - max number of three-body potentials in system
  * - ``mx2tbp``
    - *variable*
    - array dimension of three-body potential parameters
  * - ``mxptbp``
    - *variable* 
    - max number of three-body potential parameters (5)
  * - ``mxfbp`` 
    - *variable*
    - max number of four-body potentials in system
  * - ``mx2fbp`` 
    - *variable*
    - array dimension of four-body potential parameters
  * - ``mxpfbp``
    - *variable*
    - max number of four-body potential parameters (3)
  * - ``mxpfld``
    - *variable*
    - max number of external field parameters (5)
  * - ``mxstak``
    - *variable*
    - dimension of stack arrays for rolling averages
  * - ``mxnstk``
    - *variable*
    - max number of stacked variables
  * - ``mxlist``
    - *variable*
    - max number of atoms in the Verlet list on a node
  * - ``mxcell``
    - *variable* 
    - max number of link cells per node 
  * - ``mxatms`` 
    - *variable*
    - max number of local+halo atoms per node
  * - ``mxatdm`` 
    - *variable*
    - max number of local atoms per node
  * - ``mxbfdp``
    - *variable* 
    - max dimension of the transfer buffer for deport functions
  * - ``mxbfss`` 
    - *variable* 
    - max dimension of the transfer buffer for statistics functions
  * - ``mxbfxp`` 
    - *variable* 
    - max dimension of the transfer buffer for export functions
  * - ``mxbfsh`` 
    - *variable*
    - max dimension of the transfer buffer for shared units
  * - ``mxbuff``
    - *variable*  
    - max dimension of the principle transfer buffer
  * - ``zero_plus``
    - *variable* 
    - the machine representation of :math:`+0` at working precision
  * - ``half_plus`` 
    - *variable* 
    - the machine representation of :math:`+0.5\uparrow` at working precision
  * - ``half_minus`` 
    - *variable* 
    - the machine representation of :math:`+0.5\downarrow` at working precision
