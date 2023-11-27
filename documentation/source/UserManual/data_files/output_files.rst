.. _output-files:

The OUTPUT Files
================

DL_POLY_4 may produce many output files. However only OUTPUT (an
incremental summary file of the simulation), STATIS (a statistical
history file), REVCON (a restart configuration file - final) and REVIVE
(a restart statistics accumulators file - final) are mandatory. DUMP_E
(a restart electronic temperature grid file - final) is also produced if
the two-temperature model (\ :index:`TTM<Two-Temperature Model>`) 
is in use. The existence of the
remaining files is optional upon user specifications in CONTROL. Some of
these optional files are HISTORY, DEFECTS, MSDTMP, CFGMIN, RDFDAT,
USRDAT, ZDNDAT, VDFDAT, LATS_E, LATS_I, PEAK_E, PEAK_I. These
respectively contain: an incremental dump file of all atomic
coordinates, velocities and forces; an incremental dump file of atomic
coordinates of defected particles (interstitials) and sites (vacancies);
an incremental dump file of of individual atomic mean square
displacement and temperature; a dump file of all atomic coordinates of a
minimised structure; a radial distribution function (RDF) data file; the
RDF data file for the umbrella sampling (harmonic restraint); Z-density
distribution data file; velocity autocorrelation function (VAF) data
files (one file for each species); electronic temperature profile data
file; ionic temperature profile data file; electronic temperature
statistical data file; ionic temperature statistical data file.

.. _history-file:

The HISTORY File
----------------

The HISTORY file is the dump file of atomic coordinates, velocities and
forces. Its principal use is for off-line analysis. The file is written
by the subroutine ``trajectory_write``. The control variables for this
file are ``ltraj, nstraj, istraj`` and keytrj which are created
internally, based on information read from the **traj** directive in the
CONTROL file (see Section \ :ref:`control-file`). The
HISTORY file will be created only if the directive **traj** appears in
the CONTROL file.

The HISTORY file can become *very* large, especially if it is formatted.
For serious simulation work it is recommended that the file be written
to a scratch disk capable of accommodating a large data file.
Alternatively, the file may be written in netCDF format instead of in
ASCII (users must change ensure this functionality is available), which
has the additional advantage of speed.

The HISTORY has the following structure:

.. parsed-literal:: 

  **record 1**
  ``header``    a72       file header
  **record 2**
  ``keytrj``    integer   trajectory key (see Table \ `[keytrj] <#keytrj>`__) in last frame
  ``imcon``     integer   periodic boundary key (see Table :numref:`(%s)<imcon_table>`) in last frame
  ``megatm``    integer   number of atoms in simulation cell in last frame
  ``frame``     integer   number configuration frames in file
  ``records``   integer   number of records in file

For timesteps greater than ``nstraj`` the HISTORY file is appended at
intervals specified by the **traj** directive in the CONTROL file, with
the following information for each configuration:

.. parsed-literal:: 
  
  **record i**
  ``timestep``    a8        the character string “timestep”
  ``nstep``       integer   the current time-step
  ``megatm``      integer   number of atoms in simulation cell (again)
  ``keytrj``      integer   trajectory key (again)
  ``imcon``       integer   periodic boundary key (again)
  ``tstep``       real      integration timestep (ps)
  ``time``        real      elapsed simulation time (ps)
  **record ii**
  ``cell(1)``     real      x component of *a* cell vector in Å
  ``cell(2)``     real      y component of *a* cell vector in Å
  ``cell(3)``     real      z component of *a* cell vector in Å
  **record iii**
  ``cell(4)``     real      x component of *b* cell vector in Å
  ``cell(5)``     real      y component of *b* cell vector in Å
  ``cell(6)``     real      z component of *b* cell vector in Å
  **record iv**
  ``cell(7)``     real      x component of *c* cell vector in Å
  ``cell(8)``     real      y component of *c* cell vector in Å
  ``cell(9)``     real      z component of *c* cell vector in Å

This is followed by the configuration for the current timestep. i.e. for
each atom in the system the following data are included:

.. parsed-literal:: 
  
  **record a**
  ``atmnam``    a8        atomic label
  ``iatm``      integer   atom index
  ``weight``    real      atomic mass (a.m.u.)
  ``charge``    real      atomic charge (e)
  ``rsd``       real      displacement from position at *t* = 0 in Å
  **record b** 
  ``xxx``       real      x coordinate
  ``yyy``       real      y coordinate
  ``zzz``       real      z coordinate
  **record c** only for ``keytrj`` > 0
  ``vxx``       real      x component of velocity in Å/picosecond
  ``vyy``       real      y component of velocity in Å/picosecond
  ``vzz``       real      z component of velocity in Å/picosecond
  **record d** only for ``keytrj`` > 1
  ``fxx``       real      x component of force in Å\ :math:`\cdot`\ Dalton/picosecond\ :math:`^{2}`
  ``fyy``       real      y component of force in Å\ :math:`\cdot`\ Dalton/picosecond\ :math:`^{2}`
  ``fzz``       real      z component of force in Å\ :math:`\cdot`\ Dalton/picosecond\ :math:`^{2}`

Thus the data for each atom is a minimum of two records and a maximum of
4.

.. _msdtmp-file:

The MSDTMP File
---------------

The MSDTMP file is the dump file of individual atomic mean square
displacements (square roots in Å) and mean square temperature (square
roots in Kelvin). Its principal use is for off-line analysis. The file
is written by the subroutine ``msd_write``. The control variables for
this file are ``l_msd, nstmsd, istmsd`` which are created internally,
based on information read from the **msdtmp** directive in the CONTROL
file (see Section \ :ref:`control-file`). The MSDTMP file
will be created only if the directive **msdtmp** appears in the CONTROL
file.

The MSDTMP file can become *very* large, especially if it is formatted.
For serious simulation work it is recommended that the file be written
to a scratch disk capable of accommodating a large data file.

The MSDTMP has the following structure:

.. parsed-literal::

  **record 1**
  ``header``    a52       file header
  **record 2**
  ``megatm``    integer   number of atoms in simulation cell in last frame
  ``frame``     integer   number configuration frames in file
  ``records``   integer   number of records in file

For timesteps greater than ``nstmsd`` the MSDTMP file is appended at
intervals specified by the **msdtmp** directive in the CONTROL file,
with the following information for each configuration:

.. parsed-literal::

  **record i**
  ``timestep``    a8 the      character string “timestep”
  ``nstep``       integer     the current time-step
  ``megatm``      integer     number of atoms in simulation cell (again)
  ``tstep``       real        integration timestep (ps)
  ``time``        real        elapsed simulation time (ps)

This is followed by the configuration for the current timestep. i.e. for
each atom in the system the following data are included:

.. parsed-literal::

  **record a**
  ``atmnam``    a8        atomic label
  ``iatm``      integer   atom index :math:`\sqrt{\texttt{MSD}(t)}` real square root of the atomic mean square displacements (in Å)
  T(mean)       real      atomic mean temperature (in Kelvin)

.. _defects-file:

The DEFECTS File
----------------

The DEFECTS file is the dump file of atomic coordinates of defects (see
Section \ `The REFERENCE File<reference-file>`). Its principal use is
for off-line analysis. The file is written by the subroutine
``defects_write``. The control variables for this file are
``ldef, nsdef, isdef`` and ``rdef`` which are created internally, based
on information read from the **defects** directive in the CONTROL file
(see Section \ :ref:`control-file`). The DEFECTS file
will be created only if the directive **defects** appears in the CONTROL
file.

The DEFECTS file may become *very* large, especially if it is formatted.
For serious simulation work it is recommended that the file be written
to a scratch disk capable of accommodating a large data file.

The DEFECTS has the following structure:

.. parsed-literal::
  
  **record 1**
  ``header``    a72       file header
  **record 2**
  ``rdef``      real      site-interstitial cutoff (Å) in last frame
  ``frame``     integer   number configuration frames in file
  ``records``   integer   number of records in file

For timesteps greater than ``nsdef`` the DEFECTS file is appended at
intervals specified by the **defects** directive in the CONTROL file,
with the following information for each configuration:

.. parsed-literal:: 
  
  **record i**
  ``timestep``      a8        the character string “timestep”
  ``nstep``         integer   the current time-step
  ``tstep``         real      integration timestep (ps)
  ``time``          real      elapsed simulation time (ps)
  ``imcon``         integer   periodic boundary key (see Table :numref:`(%s)<imcon_table>`)
  ``rdef``          real      site-interstitial cutoff (Å)
  **record ii**
  ``defects``       a7        the  character string “defects”
  ``ndefs``         integer   the total number of defects
  ``interstitials`` a13       the character string “interstitials”
  ``ni``            integer   the total number of interstitials
  ``vacancies``     a9        the character string “vacancies”
  ``nv``            integer   the total number of vacancies
  **record iii**
  ``cell(1)``       real      x component of *a* cell vector
  ``cell(2)``       real      y component of *a* cell vector
  ``cell(3)``       real      z component of *a* cell vector
  **record iv**
  ``cell(4)``       real      x component of *b* cell vector
  ``cell(5)``       real      y component of *b* cell vector
  ``cell(6)``       real      z component of *b* cell vector
  **record v**
  ``cell(7)``       real      x component of *c* cell vector
  ``cell(8)``       real      y component of *c* cell vector
  ``cell(9)``       real      z component of *c* cell vector

This is followed by the ``ni`` interstitials for the current timestep,
as each interstitial has the following data lines:

.. parsed-literal:: 
  
  **record a**
  ``atmnam``    a10       i_atomic label from CONFIG
  ``iatm``      integer   atom index from CONFIG
  **record b**
  ``xxx``       real      x coordinate
  ``yyy``       real      y coordinate
  ``zzz``       real      z coordinate

This is followed by the ``nv`` vacancies for the current timestep, as
each vacancy has the following data lines:

.. parsed-literal:: 
  
  **record a**
  ``atmnam``  a10       v_atomic label from REFERENCE
  ``iatm``    integer   atom index from REFERENCE
  **record b**
  ``xxx``     real      x coordinate from REFERENCE
  ``yyy``     real      y coordinate from REFERENCE
  ``zzz``     real      z coordinate from REFERENCE

.. _rsddat-file:

The RSDDAT File
---------------

The RSDDAT file is the dump file of atomic coordinates of atoms that are
displaced from their original position at :math:`t~=~0` farther than a
preset cutoff. Its principal use is for off-line analysis. The file is
written by the subroutine ``rsd_write``. The control variables for this
file are ``lrsd, nsrsd, isrsd`` and ``rrsd`` which are created
internally, based on information read from the **displacements**
directive in the CONTROL file (see
Section \ :ref:`control-file`). The RSDDAT file will be
created only if the directive **defects** appears in the CONTROL file.

The RSDDAT file may become *very* large, especially if it is formatted.
For serious simulation work it is recommended that the file be written
to a scratch disk capable of accommodating a large data file.

The RSDDAT has the following structure:

.. parsed-literal:: 
  
  **record 1**
  ``header``    a72       file header
  **record 2**
  ``rdef``      real      displacement qualifying cutoff (Å) in last frame
  ``frame``     integer   number configuration frames in file
  ``records``   integer   number of records in file

For timesteps greater than ``nsrsd`` the RSDDAT file is appended at
intervals specified by the **displacements** directive in the CONTROL
file, with the following information for each configuration:

.. parsed-literal:: 
  
  **record i**
  ``timestep``        a8        the character string “timestep”
  ``nstep``           integer   the current time-step
  ``tstep``           real      integration timestep (ps)
  ``time``            real      elapsed simulation time (ps)
  ``imcon``           integer   periodic boundary key (see Table :numref:`(%s)<imcon_table>`)
  ``rrsd``            real      displacement qualifying cutoff (Å)
  **record ii**
  ``displacements``   a13       the character string “displacements”
  ``nrsd``            integer   the total number of displacements
  **record iii**
  ``cell(1)``         real      x component of *a* cell vector
  ``cell(2)``         real      y component of *a* cell vector
  ``cell(3)``         real      z component of *a* cell vector
  **record iv**
  ``cell(4)``         real      x component of *b* cell vector
  ``cell(5)``         real      y component of *b* cell vector
  ``cell(6)``         real      z component of *b* cell vector
  **record v**
  ``cell(7)``         real      x component of *c* cell vector
  ``cell(8)``         real      y component of *c* cell vector
  ``cell(9)``         real      z component of *c* cell vector

This is followed by the ``nrsd`` displacements for the current timestep,
as each atom has the following data lines:

.. parsed-literal:: 
  
  **record a**
  ``atmnam``    a10       atomic label from CONFIG
  ``iatm``      integer   atom index from CONFIG
  ``ratm``      real      atom displacement from its position at :math:`t~=~0`
  **record b**
  ``xxx``       real      x coordinate
  ``yyy``       real      y coordinate
  ``zzz``       real      z coordinate

.. _cfgminfile:

The CFGMIN File
---------------

The CFGMIN file only appears if the user has selected the programmed
minimisation option (directive **minim**\ ise (or **optim**\ ise) in the
CONTROL file). Its contents have the same format as the CONFIG file (see
Section \ :ref:`config-file`), but contains only atomic
position data and will never contain either velocity or force data (i.e.
parameter ``levcfg`` is always zero). In addition, three extra numbers
appear on the end of the second line of the file:

#. an integer indicating the number of minimisation cycles required to
   obtain the structure,

#. the configuration energy of the minimised configuration expressed in
   DL_POLY_4 units (Section `[units] <#units>`__), and

#. the configuration energy of the initial structure expressed in
   DL_POLY_4 units (Section `[units] <#units>`__).

.. _output-file:

The OUTPUT File
---------------

The job output consists of 7 sections: Header; Simulation control
specifications; Force field specification; System specification; Summary
of the initial configuration; Simulation progress; Sample of the final
configuration; Summary of statistical data; and Radial distribution
functions and Z-density profile. These sections are written by different
subroutines at various stages of a job. Creation of the OUTPUT file
*always* results from running . It is meant to be a human readable file,
destined for hardcopy output.

Header
~~~~~~

Gives the DL_POLY_4 version number, the number of processors in use, the
link-cell algorithm in use and a title for the job as given in the
header line of the input file CONTROL. This part of the file is written
from the subroutines ``dl_poly``, ``set_bounds`` and ``read_control``.

Simulation Control Specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Echoes the input from the CONTROL file. Some variables may be reset if
illegal values were specified in the CONTROL file. This part of the file
is written from the subroutine read_control.

Force Field Specification
~~~~~~~~~~~~~~~~~~~~~~~~~

Echoes the FIELD file. A warning line will be printed if the system is
not electrically neutral. This warning will appear immediately before
the non-bonded short-range potential specifications. This part of the
file is written from the subroutine ``read_field``.

System Specification
~~~~~~~~~~~~~~~~~~~~

Echoes system name, periodic boundary specification, the cell vectors
and volume, some initial estimates of long-ranged corrections the energy
and pressure (if appropriate), some concise information on topology and
degrees of freedom break-down list. This part of the file is written
from the subroutines ``scan_config``, ``check_config``, ``system_init``,
``report_topology`` and ``set_temperature``.

Summary of the Initial Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This part of the file is written from the main subroutine dl_poly\_. It
states the initial configuration of (a maximum of) 20 atoms in the
system. The configuration information given is based on the value of
``levcfg`` in the CONFIG file. If levcfg is 0 (or 1) positions (and
velocities) of the 20 atoms are listed. If ``levcfg`` is 2 forces are
also written out.

Simulation Progress
~~~~~~~~~~~~~~~~~~~

This part of the file is written by the DL_POLY_4 root segment
``dl_poly``. The header line is printed at the top of each page as:

::

   --------------------------------------------------------------------------------------------------

       step   eng_tot  temp_tot   eng_cfg   eng_src   eng_cou   eng_bnd   eng_ang   eng_dih   eng_tet
   time(ps)    eng_pv  temp_rot   vir_cfg   vir_src   vir_cou   vir_bnd   vir_ang   vir_con   vir_tet
   cpu  (s)    volume  temp_shl   eng_shl   vir_shl     alpha      beta     gamma   vir_pmf     press

   --------------------------------------------------------------------------------------------------

The labels refer to :

.. parsed-literal:: 
  
  **line 1**
  ``step``         MD step number
  ``eng_tot``     total internal energy of the system
  ``temp_tot``    system temperature (in Kelvin)
  ``eng_cfg``     configurational energy of the system
  ``eng_src``     configurational energy due to short-range potential contributions
  ``eng_cou``     configurational energy due to electrostatic :index:`potential<potential;electrostatics>`
  ``eng_bnd``     configurational energy due to chemical bond :index:`potentials<potential;bond>`
  ``eng_ang``     configurational energy due to :index:`valence<potential;valence angle>` angle and :index:`three-body<potential;three-body>` potentials
  ``eng_dih``     configurational energy due to :index:`dihedral<potential;dihedral>` inversion and :index:`four-body<potential;four-body>` potentials
  ``eng_tet``     configurational energy due to :index:`tethering<potential;tether>` potentials
  **line 2**
  ``time(ps)``    elapsed simulation time (in pico-seconds) since the beginning of the job
  ``eng_pv``      enthalpy of system
  ``temp_rot``    rotational temperature (in Kelvin)
  ``vir_cfg``     total configurational contribution to the virial
  ``vir_src``     short range potential contribution to the virial
  ``vir_cou``     electrostatic :index:`potential<potential;electrostatics>` contribution to the virial
  ``vir_bnd``     chemical bond contribution to the virial
  ``vir_ang``     angular and :index:`three-body<potential;three-body>` potentials contribution to the virial
  ``vir_con``     constraint bond contribution to the virial
  ``vir_tet``     tethering :index:`potential<potential;tether>` contribution to the virial
  **line 3**
  ``cpu (s)``     elapsed cpu time (in seconds) since the beginning of the job
  ``volume``      system volume (in Å\ :math:`^{3}`)
  ``temp_shl``    core-shell temperature (in Kelvin)
  ``eng_shl``     configurational energy due to core-shell potentials
  ``vir_shl``     core-shell potential contribution to the virial
  ``alpha``       angle between *b* and *c* cell vectors (in degrees)
  ``beta``        angle between *c* and *a* cell vectors (in degrees)
  ``gamma``       angle between *a* and *b* cell vectors (in  degrees)
  ``vir_pmf``     :index:`PMF<constraints;PMF>` constraint contribution to the virial
  ``press``       pressure (in kilo-atmospheres)

**Note:** The total internal energy of the system (variable
``tot_energy``) includes all contributions to the energy (including
system extensions due to thermostats etc.). It is nominally the
*conserved variable* of the system, and is not to be confused with
conventional system energy, which is a sum of the kinetic and
configuration energies.

The interval for printing out these data is determined by the directive
**print** in the CONTROL file. At each time-step that printout is
requested the instantaneous values of the above statistical variables
are given in the appropriate columns. Immediately below these three
lines of output the rolling averages of the same variables are also
given. The maximum number of time-steps used to calculate the rolling
averages is controlled by the directive **stack** in file CONTROL (see
above) and listed as parameter ``mxstak`` in the ``setup_module`` file
(see Section :ref:`file-structure`). The default
value is ``mxstak`` :math:`=~100`.

Energy Units
~~~~~~~~~~~~

The energy unit for the energy and virial data appearing in the OUTPUT
is defined by the **units** directive appearing in the FIELD file.
System energies are therefore read in **units** per MD cell.

Pressure Units
~~~~~~~~~~~~~~

.. index:: single: units;pressure

The unit of pressure is katms, irrespective of what energy unit is
chosen.

Two-Temperature Model
~~~~~~~~~~~~~~~~~~~~~

.. index:: single: Two-Temperature Model

If the two-temperature model is in use, information about the timestep
sizes used for electronic thermal diffusivity is written immediately
prior to each report of statistical variables at each molecular dynamics
timestep for which printout is requested. The optimum diffusive timestep
size is given in pico-seconds, along with the chosen value and the
corresponding number of divisions of the MD timestep. If dynamic
calculation of the average atomic density in active cells is requested,
this value is included along with the number of active ionic temperature
cells. Reports are also given when energy deposition starts and
finishes.

Sample of Final Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The positions, velocities and forces of the 20 atoms used for the sample
of the initial configuration (see above) are given. This is written by
the main subroutine ``dl_poly``.

Summary of Statistical Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This portion of the OUTPUT file is written from the subroutine
``statistics_result``. The number of time-steps used in the collection
of statistics is given. Then the averages over the production portion of
the run are given for the variables described in the previous section.
The root mean square variation in these variables follow on the next two
lines. The :index:`energy<units;DL_POLY>` and :index:`pressure<units;pressure>` 
units are as for the preceding section.

Also provided in this section are estimates of the diffusion coefficient
and the mean square displacement for the different atomic species in the
simulation. These are determined from a single time origin and are
therefore approximate. Accurate determinations of the diffusion
coefficients can be obtained using the ``msd`` utility program, which
processes the HISTORY file (see User Manual).

If an NPT (N:math:`\mat{\sigma}`\ T) simulation is performed the OUTPUT
file also provides the mean pressure (and stress tensor in pressure
units as density) and mean simulation cell vectors. In case when
extended N\ :math:`\underline{\underline{\mathbf{\sigma}}}`\ T ensembles are used then further mean
:math:`(x,y)` plain area and mean surface tension are also displayed in
the OUTPUT file.

Radial Distribution Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If both calculation and printing of radial distribution functions have
been requested (by selecting directives **rdf** and print rdf in the
CONTROL file) radial distribution functions are printed out. This is
written from the subroutine rdf_compute. First the number of time-steps
used for the collection of the histograms is stated.

For each function a header line states the atom types (‘a’ and ‘b’)
represented by the function. Then :math:`r, g(r)` and :math:`n(r)` are
given in tabular form. :math:`n(r)` is the average number of atoms of
type ‘b’ within a sphere of radius :math:`r` around an atom of type ‘a’.
Note that a readable version of these data is provided by the RDFDAT
file (below).

Umbrella Sampling Restraint RDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If an umbrella sampling harmonic restraint is defined in the FIELD file
(by selecting the **ushr** external field sectione) the RDF of the two
restraint objects/fragments is printed out. This is written from the
subroutine ``usr_compute`` in ``rdf_compute``. Note that a readable
version of these data is provided by the USRDAT file (below).

Z-density Profile
~~~~~~~~~~~~~~~~~

If both calculation and printing of Z-density profiles have been
requested (by selecting directives **zden** and **print zden** in the
CONTROL file Z-density profiles are printed out as the last part of the
OUTPUT file. This is written by the subroutine z_density_compute. First
the number of time-steps used for the collection of the histograms is
stated. Then each function is given in turn. For each function a header
line states the atom type represented by the function. Then
:math:`z,~\rho(z)` and :math:`n(z)` are given in tabular form. Output is
given from :math:`Z = [-L/2,L/2]` where L is the length of the MD cell
in the Z direction and :math:`\rho(z)` is the mean number density.
:math:`n(z)` is the running integral from :math:`-L/2` to :math:`z` of
:math:`({\rm xy~cell~area}) \times \rho(s)~ds`. Note that a readable
version of these data is provided by the ZDNDAT file (below).

Velocity Autocorrelation Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If both calculation and printing of velocity autocorrelation functions
have been requested (by selecting directives **vaf** and **print vaf**
in the CONTROL file the velocity autocorrelation function for the system
(either time-averaged or the last complete sample) is printed out as the
last part of the OUTPUT file. This is written by the subroutine
``vaf_compute``. First the details of the calculations are stated:
either the number of samples used to give a time-averaged profile or the
number of the last completed sample with its starting time. The absolute
value of the velocity autocorrelation function for the system at
:math:`t=0`, :math:`C(0)`, is then stated. Then :math:`t` and
:math:`Z(t)` are given in tabular form. :math:`Z(t)=C(t)/C(0)` is the
value of the velocity autocorrelation function,
:math:`C(t)=\langle \underline{v}_{i}(0) \cdot  \underline{v}_{i}(t) \rangle`,
scaled by :math:`C(0) \equiv 3k_B T/m`. Note that a readable version of
these data for individual species is provided by the VAFDAT files
(below).

The HEATFLUX File
-----------------

The HEATFLUX file contains data relevant to the calculation of heat-flux
via a Green-Kubo mothod via an external convolution, the information is
written as:

::

   STEP  STPTMP   VOLUME   HEAT_FLUX

.. _pp-cont-file:

The PP_CONT File
----------------

This file contains the contributions of each particle to energies,
forces and stresses in a format similar to to the CONFIG file, but with
ID replaced with energy, and velocities/forces with the stress 6-vector.

::

     TAG  ATMNAM  KIN_E  MASS ENERGY
     STR_XX  STR_YX  STR_ZX
     STR_XY  STR_YY  STR_ZY
     STR_XZ  STR_YZ  STR_ZZ

.. _revcon-file:

The REVCON File
---------------

This file is formatted and written by the subroutine ``revive``. REVCON
is the restart configuration file. The file is written every ``ndump``
time steps in case of a system crash during execution and at the
termination of the job. A successful run of DL_POLY_4 will always
produce a REVCON file, but a failed job may not produce the file if an
insufficient number of timesteps have elapsed. ndump is controlled by
the directive **dump** in file CONTROL (see above) and listed as
parameter ``ndump`` in the setup_module file (see
Section :ref:`file-structure`). The default value is
``ndump`` :math:`=1000`. REVCON is identical in format to the CONFIG
input file (see Section \ :ref:`config-file`). REVCON
should be renamed CONFIG to continue a simulation from one job to the
next. This is done for you by the copy macro supplied in the execute
directory of .

.. _revive-file:

The REVIVE File
---------------

This file is unformatted and written by the subroutine system_revive. It
contains the accumulated statistical data. It is updated whenever the
file REVCON is updated (see previous section). REVIVE should be renamed
REVOLD to continue a simulation from one job to the next. This is done
by the copy macro supplied in the *execute* directory of . In addition,
to continue a simulation from a previous job the restart keyword must be
included in the CONTROL file.

The format of the REVIVE file is identical to the REVOLD file described
in Section \ :ref:`revold-file`.

The DUMP_E File
---------------

This file is formatted and written by the subroutine
``ttm_system_revive`` every ``ndump`` time steps. It contains the
electronic temperatures of all coarse-grained electronic temperature
(CET) cells and can be used to restart a simulation using the
two-temperature :index:`model<Two-Temperature Model>` without renaming 
the file.

The format of the DUMP_E is described in
Section \ :ref:`dumpe-file`.

.. _rdf-file:

The RDFDAT File
---------------

This is a formatted file containing *Radial Distribution Function* (RDF)
data. Its contents are as follows:

.. parsed-literal::
  
  **record 1**
  ``cfgname``   a72       configuration name
  **record 2**
  ``ntprdf``    integer   number of different RDF pairs tabulated in file
  ``mxgrdf``    integer   number of grid points for each RDF pair

There follow the data for each individual RDF, i.e. ``ntprdf`` times.
The data supplied are as follows:

.. parsed-literal::
  
  **first record**
  ``atname 1``  a8      first atom name
  ``atname 2``  a8      second atom name
  **following records** (*mxgrdf* records)
  ``radius``    real    interatomic distance (Å)
  ``g(r)``      real    RDF at given radius
  ``n(r)``      real    RDF at given radius

**Note 1.** The RDFDAT file is optional and appears when the **print
rdf** option is specified in the CONTROL file.

**Note 2.** Along with the RDFDAT file, two other files will be created
whenever the **print ana**\ lysis directive is invoked: VDWPMF & VDWTAB,
both containing the data for potentials of mean force and the
corresponding virials calculated based on the obtained RDF:s, i.e. PMF
:math:`\sim -\ln({\rm RDF})` (in the energy units specified in the FIELD
file). These files have a simple three column format, the same as that
used for \*PMF files in the case of bonded units, see
Section \ :ref:`IPDF-analysis`. The purpose of these
files is to provide the user with means of setting up a PMF-based
force-field, for example in the case of initial coarse-graining of an
atomistic system. In particular, one can convert the VDWTAB file into a
correctly formatted TABLE file (Section :ref:`table-file`)
by using the utility called ``pmf2tab.f`` (subject to compilation; found
in DL_POLY_4 directory ``utility``) as follows,

``[user@host]$ pmf2tab.exe < VDWTAB``

see Section \ :ref:`cg-intro` for completeness.

.. _usr-file:

The USRDAT File
---------------

.. parsed-literal:: 
  
  **record 1**
  ``# title``   a100      file header title
  **record 2**
  ``# header``  a100      file information header
  **record 3**
  ``# info``    a30       information to follow string
  **record 3**
  ``bins``      integer   number of bins
  ``cutoff``    real      cutoff in Å
  ``frames``    integer   number of sampled configurations
  ``volume``    real      average cell volue Å\ :math:`^{3}`
  **record 4**
  ``#``         a1        a hash (#) symbol
  **following records** (*mxgusr* records)
  ``radius``    real      interatomic distance (Å)
  ``g(r)``      real      RDF at given radius

.. _zdn-file:

The ZDNDAT File
---------------

This is a formatted file containing the Z-density data. Its contents are
as follows:

.. parsed-literal:: 
  
  **record 1**
  ``cfgname``   a72       configuration name
  **record 2**
  ``ntpatm``    integer   number of unique atom types profiled in file
  ``mxgrdf``    integer   number of grid points in the Z-density function 

There follow the data for each individual Z-density function, i.e.
``ntpatm`` times. The data supplied are as follows:

.. parsed-literal::
  
  **first record**
  ``atname``    a8    unique atom name
  **following records** (*mxgrdf* records)
  ``z``         real  distance in z direction (Å)
  :math:`\rho(z)` real Z-density at given height ``z``

**Note** the ZDNDAT file is optional and appears when the **print rdf**
option is specified in the CONTROL file.

.. _vaf-files:

The VAFDAT Files
----------------

These are formatted files containing *Velocity Autocorrelation Function*
(VAF) data. An individual file is created for each atomic species, i.e.
``VAFDAT``\ \_\ *atname*. Their contents are as follows:

.. parsed-literal:: 
  
  **record**
  ``cfgname``   a72   configuration name

There follow the data for the VAF, either a single time-averaged profile
or successive profiles separated by two blank lines. The data supplied
are as follows:

.. parsed-literal:: 
  
  **first record**
  ``atname``    a8        atom name
  ``binvaf``    integer   number of data points in VAF profile, *excluding* :math:`t=0`
  ``vaforigin`` real      absolute value of VAF at :math:`t=0`  (:math:`C(0) \equiv 3k_B T/m`)
  ``vaftime0``  real      simulation time (ps) at beginning of (last) VAF profile (:math:`t=0`)
  **following records** (*binvaf*\ +1 records)
  ``t``         real      time (ps)
  ``Z(t)``      real      scaled velocity autocorrelation function (:math:`C(t)/C(0)`) at given time :math:`t`

**Note** the VAFDAT files are optional and appear when the **print vaf**
option is specified in the CONTROL file.

.. _bonded-files:

The INTDAT, INTPMF & INTTAB Files
---------------------------------

These files, where INT is referring to INTra-molecular interactions and
VDW(RDF derived inter-molecular), have very similar formatting rules
with some examples shown in
Section \ :ref:`IPDF-analysis`. Refer to
Section \ :ref:`IPDF-analysis` for their meaning and
usage in coarse grained model systems.

.. parsed-literal:: 
  
  **record 1**
  ``# title``     a100        file header title
  **record 2**
  ``# header``    a100        file information header
  **record 3**
  ``# info``      a30         information to follow string
  ``bins``        integer     number of bins for all PDFs
  ``cutoff``      real        cutoff in Å for bonds and RDFs or degrees for angular intramolecular interactions
  ``frames``      integer     number of sampled configurations
  ``types``       integer     number of unique types of these interactions
  **record 4**
  ``#``           a1          a hash (#) symbol
  **record 5**
  ``# info 1``    a100        information to follow string
  **record 6**
  ``#``           a1          a hash (#) symbol

The subsequent records define each PDF potential in turn, in the order
indicated by the specification in the FIELD file. Each potential is
defined by a header record and a set of data records with the
potential-like and force-like tables.

.. parsed-literal:: 
  
  **empty record:**
  **id record:**
  ``# info``    a25     information to follow string
  ``atom 1``    a8      first atom type
  ``atom 2``    a8      second atom type
  ``atom 3``    a8      third atom type - only available in ANG\* files
  ``atom 4``    a8      forth atom type - only available in DIH\* & INV\* files
  ``index``     integer unique index of PDF in file
  ``instances`` integer instances of this unique type of PDF
  **interaction data records 1–bins:**
  ``abscissa``  real    consecutive value over the full cutoff/range in Å for BNDTAB & VDWTAB and degrees for ANGTAB, DIHTAB & INVTAB
  ``potential`` real    potential at the abscissa grid point in **units** as specified in FIELD
  ``force``     real    complementary force (virial for BNDTAB & VDWTAB) value

.. _statis-file:

The STATIS File
---------------

The file is formatted, with integers as “i10” and reals as “e14.6”. It
is written by the subroutine statistics_collect. It consists of two
header records followed by many data records of statistical data.

.. parsed-literal:: 
  
  **record 1**
  ``cfgname``   a72   configuration name
  **record 2**
  ``string``    a8    energy units

Data records
~~~~~~~~~~~~

Subsequent lines contain the instantaneous values of statistical
variables dumped from the array stpval. A specified number of entries
of ``stpval`` are written in the format “(1p,5e14.6)”. The number of
array elements required (determined by the parameter ``mxnstk`` in the
``setup_module`` file) is

  .. math::

     \begin{aligned}
     \texttt{mxnstk} \ge ~& 28 + 9~(\rm stress~tensor~elements) ~+ \nonumber \\
      & \texttt{ntpatm}~(\rm number~of~unique~atomic~sites) ~+ \nonumber \\
      & 10~(\rm if~constant~pressure~simulation~requested) ~+ \nonumber \\
      & 2~(\rm if~iso~>~0~requested) + 2~(\rm if~iso~>~1~requested) ~+ \nonumber \\
      & 2*mxatdm~(\rm if~msdtmp~option~is~used) \nonumber\end{aligned}

The STATIS file is appended at intervals determined by the stats
directive in the CONTROL file. The energy unit is as specified in the
FIELD file with the **units** directive, and are compatible with the
data appearing in the OUTPUT file. The contents of the appended
information of calculated *instantaneous* observables is:

.. parsed-literal:: 
  
  **record i**
  ``nstep``         integer   current MD time-step
  ``time``          real      elapsed simulation time
  ``nument``        integer   number of array elements to follow
  **record ii** ``stpval``\ (1) – ``stpval``\ (5)
  ``engcns``        real      total extended system energy, :math:`E^{x}_{tot}=(E_{kin}+E_{rot})+E_{conf}+E_{consv}` (i.e. including the conserved quantity, :math:`E_{consv}`)
  ``temp``          real      system temperature, :math:`2\frac{E_{kin}+E_{rot}}{f k_{B}}`
  ``engcfg``        real      configurational energy, :math:`E_{conf}`
  ``engsrc``        real      short range potential energy
  ``engcpe``        real      electrostatic energy
  **record iii** ``stpval``\ (6) – ``stpval``\ (10)
  ``engbnd``        real      chemical bond energy
  ``engang``        real      valence angle and 3-body potential energy
  ``engdih``        real      dihedral, inversion, and 4-body potential energy
  ``engtet``        real      tethering energy
  ``enthal``        real      enthalpy (:math:`E^{x}_{tot} + {\cal P} \cdot V`) for NVE/T/E\ :math:`_{kin}` ensembles
  enthalpy (:math:`E^{x}_{tot} + P \cdot {\cal V}`) for NP/\ :math:`\sigma`\ T or NP\ :math:`_{n}`\ A/\ :math:`\gamma` ensembles
  **record iv** ``stpval``\ (11) – ``stpval``\ (15)
  ``tmprot``        real      rotational temperature, :math:`E_{rot}`
  ``vir``           real      total virial
  ``virsrc``        real      short-range virial
  ``vircpe``        real      electrostatic virial
  ``virbnd``        real      bond virial
  **record v** ``stpval``\ (16) – ``stpval``\ (20)
  ``virang``        real      valence angle and 3-body virial
  ``vircon``        real      constraint bond virial
  ``virtet``        real      tethering virial
  ``volume``        real      volume, :math:`{\cal V}`
  ``tmpshl``        real      core-shell temperature
  **record vi** ``stpval``\ (21) – ``stpval``\ (25)
  ``engshl``        real      core-shell potential energy
  ``virshl``        real      core-shell virial
  ``alpha``         real      MD cell angle :math:`\alpha`
  ``beta``          real      MD cell angle :math:`\beta`
  ``gamma``         real      MD cell angle :math:`\gamma`
  **record vii** ``stpval``\ (26), ``stpval``\ (27), ``stpval``\ (0)
  ``virpmf``        real      PMF constraint virial
  ``press``         real      pressure, :math:`{\cal P}`
  ``consv``         real      extended DoF energy, :math:`E_{consv}`
  **the next 9 entries for the stress tensor in pressure units**
  ``stress(1)``     real      xx component of stress tensor
  ``stress(2)``     real      xy component of stress tensor
  ``stress(3)``     real      xz component of stress tensor
  ``stress(4)``     real      yx component of stress tensor
  ``...``           real      ...
  ``stress(9)``     real      zz component of stress tensor
  **the next ``ntpatm`` entries**
  ``amsd(1)``       real      mean squared displacement of first atom types
  ``amsd(2)``       real      mean squared displacement of second atom types
  ``...`` ... ...
  ``amsd(ntpatm)``  real      mean squared displacement of last atom types
  **the next 10 entries - if a NPT or N\ :math:`\mat{\sigma}`\ T simulation is undertaken**
  ``cell(1)``       real      x component of *a* cell vector
  ``cell(2)``       real      y component of *a* cell vector
  ``cell(3)``       real      z component of *a* cell vector
  ``cell(4)``       real      x component of *b* cell vector
  ``...``           real      ...
  ``cell(9)``       real      z component of *c* cell vector
  ``stpipv``        real      pressure, :math:`{\cal P} \cdot {\cal V}`
  **the next 2 entries - if NP\ :math:`_{n}`\ AT simulation is undertaken with iso > 0**
  ``h_z``           real      MD cell height :math:`h_{z}` to normal surface :math:`{\cal A}\perp{z}`
  ``A\perpz``       real      MD cell normal surface :math:`{\cal A}\perp{z}={\cal V}/h_{z}`
  **the next 2 entries - if a N\ :math:`\gamma_{n}`\ AT simulation is undertaken with iso > 1**
  ``gamma_x``       real      surface tension :math:`\gamma_{n_{x}}` on normal surface :math:`{\cal A}\perp{z}`
  ``gamma_y``       real      surface tension :math:`\gamma_{n_{y}}` on normal surface :math:`{\cal A}\perp{z}`



.. _latsei-files:

The LATS_E and LATS_I Files
---------------------------

These are formatted files containing electronic (LATS_E) and ionic
(LATS_I)temperatures at user-requested intervals along the y-direction
in the centre of the system’s xz-plane from two-temperature 
:index:`model<Two-Temperature Model>` calculations.

Each line in these files consists of a series of electronic or ionic
temperatures along the y-direction – ``-eltsys(2)/2`` :math:`\le y \le`
``+eltsys(2)/2`` and ``-ntsys(2)/2`` :math:`\le y \le` ``+ntsys(2)/2``
at :math:`x=z=0` – corresponding to a requested timestep. The number of
values in each line will depend on the number of electronic or ionic
temperature cells requested by the user.

.. _peakei-files:

The PEAK_E and PEAK_I Files
---------------------------

These are formatted files containing statistics from :index:`two-temperature<Two-Temperature Model>`
model calculations at user-requested intervals. Each line in these files
corresponds to a requested time step and the data is based upon active
coarse-grained electronic (CET) and ionic (CIT) temperature grid cells.

In the PEAK_E file, the data are formatted as follows:

.. parsed-literal:: 
  
  **record i**
  ``nstep``       integer current MD time-step
  ``time``        real    elapsed simulation time
  ``eltemp_min``  real    minimum value of electronic temperature in system (K)
  ``eltemp_max``  real    maximum value of electronic temperature in system (K)
  ``eltemp_mean`` real    mean value of electronic temperature in system (K)
  ``eltemp_sum``  real    sum of electronic temperatures in system (K)
  ``Ue``          real    total electronic energy in system (eV)

The PEAK_I file is formatted in a similar fashion, as follows:

.. parsed-literal:: 
  
  **record i**
  ``nstep``         integer   current MD time-step
  ``time``          real      elapsed simulation time
  ``tempion_min``   real      minimum value of ionic temperature in system (K)
  ``tempion_max``   real      maximum value of ionic temperature in system (K)
  ``tempion_mean``  real      mean value of ionic temperature in system (K)
  ``tempion_sum``   real      sum of ionic temperatures in system (K)

.. _evbpop_sec:

The POPEVB Files
----------------

This is an unformatted file to print the weight of each chemical state
:math:`|\Psi^{(k)}_{\text{EVB}}\big|^{2}` in the total EVB state, as
described in section `[sec:evb] <#sec:evb>`__. Values are printed at
each time step only after equilibration.
The structure of the printed data is as follows:
Time (ps) :math:`\,\,\,\,\,\,\,\,\,\,\,\,`
:math:`|\Psi^{(1)}_{\text{EVB}}\big|^{2}` :math:`\,\,\,\,\,\,`
:math:`|\Psi^{(2)}_{\text{EVB}}\big|^{2}` :math:`\,\,\,\,\,\,`
:math:`|\Psi^{(3)}_{\text{EVB}}\big|^{2}` :math:`\cdots\cdots`
:math:`|\Psi^{(N_F)}_{\text{EVB}}\big|^{2}`
where :math:`N_F` is the number of force-fields coupled via the EVB
simulation.

The ICOORD, CCOORD and ADFDAT files
-----------------------------------

ICOORD and CCOORD are output files that log coordination number data for
pairs of atomic species specified by the user. To perform this analysis
and output these files the user must enter the keyword
**coord_calculate** (see section
:ref:`control_options`) into the CONTROL file and
**crd** (see section :ref:`crd_sec`) into the FIELD file.

ICOORD is a dump file that can contain 2 types of data. The top of the
file contains the initial coordination of each atom and the exact atoms
it is coordinated to. There is an option to write this data at set
intervals (the writing step interval) or just at the initial step. The
bottom of the file provides the coordination distribution statistics for
each atom after each writing step interval. The coordination
distribution for the [atom list] - [atom list] pairs will also be
displayed here.

CCOORD is a coordination displacement file that dumps the positions of
all atoms that are considered to both change their initial local atomic
coordination, and move more than a set distance from their initial
position, at set intervals. This procedure is described in reference
:cite:`Diver2020`.

ADFDAT is statistics file containing the angular distributions for the
atom pairs specified.
