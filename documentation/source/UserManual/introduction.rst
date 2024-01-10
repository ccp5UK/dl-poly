.. _introduction_ch:

Introduction 
============

The DL_POLY Package
~~~~~~~~~~~~~~~~~~~

DL_POLY :cite:`smith-96a` is a package of subroutines,
programs and data files, designed to facilitate molecular dynamics
simulations of macromolecules, polymers, ionic systems and solutions on
a distributed memory parallel computer. It is available in two forms: DL_POLY_Classic
(written by Bill Smith & Tim Forester, `<http://www.ccp5.ac.uk/DL\_POLY\_CLASSIC/>`_) and DL_POLY_4(written by Ilian Todorov &
Bill Smith) :cite:`todorov-04a,todorov-06a`. Both versions
were originally written on behalf of :index:`CCP5`, the UK's Collaborative
Computational Project on Molecular Simulation, which has been in
existence since 1980 :cite:`smith-87a` (`<http://www.ccp5.ac.uk/DL\_POLY/>`_).

The two forms of DL_POLY differ primarily in their method of exploiting
parallelism. DL_POLY_Classic uses a Replicated Data (RD) strategy
:cite:`smith-91a,smith-93a,smith-94a,smith-94b` which works
well simulations of up to 30,000 atoms on up to 100 processors. DL_POLY_4 is based
on the Domain Decomposition (DD) strategy
:cite:`todorov-04a,todorov-06a,pinches-91a,rapaport-91b,smith-91a,smith-93a`,
and is best suited for large molecular simulations from :math:`10^{3}`
to :math:`10^{9}` atoms on large processor counts. The two packages are
reasonably compatible, so that it is possible to scale up from a DL_POLY_Classic to a DL_POLY_4
simulation with little effort. It should be apparent from these comments
that DL_POLY_4 is not intended as a replacement for DL_POLY_Classic.

Users are reminded that we are interested in hearing what other features
could be usefully incorporated. We obviously have ideas of our own and
CCP5 strongly influences developments, but other input would be welcome
nevertheless. We also request that our users respect the integrity of DL_POLY_4
source and not pass it on to third parties. We require that all users of
the package register with us, not least because we need to keep everyone
abreast of new developments and discovered bugs. We have developed
various forms of :index:`licence`, which we hope will ward off litigation (from
both sides), without denying access to genuine scientific users.

Further information on the DL_POLY packages may be obtained from the
DL_POLY project website -
`<http://www.ccp5.ac.uk/DL\_POLY/>`_.

Functionality
~~~~~~~~~~~~~

The following is a list of the features DL_POLY_4 supports.

Molecular Systems
-----------------

DL_POLY_4 will simulate the following molecular species:

-  Simple atomic systems and mixtures, e.g. Ne, Ar, Kr, etc.

-  Simple unpolarisable point ions, e.g. NaCl, KCl, etc.

-  Polarisable point ions and molecules, e.g. MgO, H\ :math:`_{2}`\ O,
   etc.

-  Simple :index:`rigid molecules <single: rigid body>` e.g. CCl\ :math:`_{4}`, SF\ :math:`_{6}`,
   Benzene, etc.

-  :index:`Rigid molecular<single: rigid body>` ions with point charges e.g. KNO\ :math:`_{3}`,
   (NH\ :math:`_{4}`)\ :math:`_{2}`\ SO\ :math:`_{4}`, etc.

-  Polymers with :index:`rigid bonds<single: constraints; bonds>`, e.g. C\ :math:`_{n}`\ H\ :math:`_{2n+2}`

-  Polymers with flexible and :index:`rigid bonds<single: constraints; bonds>` and point charges, e.g.
   proteins, macromolecules etc.

-  Silicate glasses and zeolites

-  Simple metals and metal alloys, e.g. Al, Ni, Cu,
   Cu\ :math:`_{3}`\ Au, etc.

-  Covalent systems as hydro-carbons and transition elements, e.g. C,
   Si, Ge, SiC, SiGe, ets.


.. index:: single: force field

Force Field
-----------

The DL_POLY_4 :index:`force field<single: force field; DL_POLY>` includes the following features:

#. All common forms of :index:`non-bonded<single: potential; non-bonded>` atom-atom (van der Waals) potentials

#. Atom-atom (and site-site) :index:`coulombic<single:potential; electrostatics>` potentials

#. :index:`Metal-metal<single: potential; metal>` (local density dependent) potentials :cite:`baskes-84a,baskes-86a,finnis-84a,sutton-90a,sutton-91a,todd-93a`

#. :index:`Tersoff<single: potential; Tersoff>` (local density dependent) potentials (for hydro-carbons) :cite:`tersoff-89a`

#. :index:`Three-body<single: potential; three-body>` :index:`valence angle<single: potential; valence angle>` and :index:`hydrogen bond<single: potential; bond>` potentials

#. :index:`Four-body<single: potential; four-body>` :index:`inversion<single: potential; inversion>` potentials

#. Ion core-shell :index:`polarasation<single: polarisation; shell model>`

#. :index:`Tether<single: potential; tether>` potentials

#. :index:`Chemical<potential; chemical bond>` bond potentials

#. :index:`Valence<potential; valence angle>` angle potentials

#. :index:`Dihedral<potential; dihedral>` angle (and :index:`improper<potential; improper dihedral>` dihedral angle) potentials

#. :index:`Inversion<potential; inversion>` angle potentials

#. :index:`External<potential; external field>` field potentials.

The parameters describing these potentials may be obtained, for example,
from the :index:`GROMOS<pair: force field; GROMOS>` :cite:`gunsteren-87a`, :index:`Dreiding<force field; Dreiding>` 
:cite:`mayo-90a` or :index:`AMBER<pair: force field; AMBER>` :cite:`weiner-86a`
forcefield, which share functional forms. It is relatively easy to adapt
to user specific :index:`force field`\ s.


.. index:: single: boundary conditions

Boundary Conditions
-------------------

DL_POLY_4 will accommodate the following boundary conditions:

#. None, e.g. isolated molecules in vacuo

#. Cubic periodic boundaries

#. Orthorhombic periodic boundaries

#. Parallelepiped periodic boundaries

#. Slab (x,y periodic, z non-periodic).

These are described in detail in 
:ref:`Appendix B<boundary-conditions>`. Note that periodic
boundary conditions (PBC) :math:`1` and :math:`5` above require careful
consideration to enable efficient load balancing on a parallel computer.


.. index:: single: Java GUI

Java Graphical User Interface
-----------------------------

.. |reg|   unicode:: U+000AE 

The DL_POLY_4 Graphical User Interface (GUI) is the same one that also comes with
DL_POLY_Classic, which is written in the Javaprogramming\ |reg| language from Sun\ |reg| Microsystems.
A major advantage of this is the free availability of the Java
programming environment from Sun\ |reg|, and also its portability across
platforms. The compiled GUI may be run without recompiling on any Java\ |reg|
supported machine. The GUI is an integral component of the DL_POLY
suites and is available on the same terms (see the GUI manual :cite:`smith-gui`).

Algorithms
----------

.. index:: 
   single: parallelisation
   single: parallelisation; Domain Decomposition

Parallel Algorithms
+++++++++++++++++++

DL_POLY_4 exclusively employs the Domain Decomposition parallelisation strategy
:cite:`pinches-91a,rapaport-91b,smith-91a,smith-93a` (see
Section: :ref:`parallelisation`).

Molecular Dynamics Algorithms
+++++++++++++++++++++++++++++

.. index:: 
   single: algorithm 
   single: algorithm; Verlet 
   single: ensemble 
   single: thermostat 
   single: barostat
   single: algorithm; RATTLE 
   single: algorithm; SHAKE
   single: constraints; bond
   single: rigid body
   single: algorithm; NOSQUISH

DL_POLY_4 offers a selection of MD integration algorithms based on Velocity Verlet
(VV) :cite:`allen-89a`. These generate NVE,
NVE\ :math:`_{kin}`, NVT, NPT and NT ensembles with a selection of
thermostats and barostats. Parallel versions of the RATTLE
:cite:`andersen-83a` and SHAKE :cite:`smith-94b`
algorithms are used for solving bond constraints. The rotational motion
of rigid bodies (RBs) is handled with the “NOSQUISH” algorithm of Miller
et al :cite:`miller-02a`.

The following MD :index:`algorithms<algorithm>` are available:

#. Constant E :index:`algorithm<ensemble;NVE>`

#. Evans constant E\ :math:`_{kin}` :index:`algorithm<ensemble;Evans NVT>` :cite:`evans-84a`

#. Langevin constant T :indeX:`algorithm<ensemble;Langevin NVT>` :cite:`adelman-76a`

#. Andersen constant T :index:`algorithm<ensemble;Andersen NVT>` :cite:`andersen-79a`

#. Berendsen constant T :index:`algorithm<ensemble;Berendsen NVT>` :cite:`berendsen-84a`

#. Nosé-Hoover constant T :index:`algorithm<ensemble;Nosé-Hoover NVT>` :cite:`hoover-85a`

#. Langevin constant T,P :index:`algorithm<ensemble;Langevin NPT>` :cite:`quigley-04a`

#. Berendsen constant T,P :index:`algorithm<ensemble;Berendsen NPT>` :cite:`berendsen-84a`

#. Nosé-Hoover constant T,P :index:`algorithm<ensemble;Nosé-Hoover NPT>` :cite:`hoover-85a`

#. Martyna, Tuckerman and Klein (MTK) constant T,P :index:`algorithm<ensemble;Martyna-Tuckerman-Klein NPT>`  :cite:`martyna-96a`

.. |sigma| unicode:: U+03A3

#. Langevin constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Langevin N|sigma|T>` :cite:`quigley-04a`

#. Berendsen constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Berendsen N$\sigma$T>` :cite:`berendsen-84a`

#. Nosé-Hoover constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Nosé-Hoover N$\sigma$T>` :cite:`hoover-85a`

#. Martyna, Tuckerman and Klein (MTK) constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Martyna-Tuckerman-Klein N$\sigma$T>` :cite:`martyna-96a`.



.. index:: dlpoly2

DL_POLY_Classic features incompatible or unavalable in DL_POLY_4
----------------------------------------------------------------

-  Force field

   -  Rigid bodies connected with constraint links are not available

   -  Shell models specification is solely determined by the presence of
      mass on the shells

   -  Dihedral potentials with more than three original parameters (see
      OPLS) have two artificially added parameters, defining the 1-4
      electrostatic and van der Waals scaling factors, which must be
      placed at fourth and fifth position respectively, extending the
      original parameter list split by them

-  Boundary conditions

   -  Truncated octahedral periodic boundaries (``imcon`` :math:`= 4`) are not
      available

   -  Rhombic dodecahedral periodic boundaries (``imcon`` :math:`= 5`) are not
      available

   -  Hexagonal prism periodic boundaries (``imcon`` :math:`= 7`) are not available

-  Electrostatics

   -  Standard Ewald Summation is not available, but is substituted by
      Smoothed Particle Mesh Ewald (SPME) summation

   -  Hautman-Klein Ewald Summation for 3D non-periodic but 2D periodic
      systems is not available

-  Non-standard functionality

   -  Temperature Accelerated Dynamics

   -  Hyperdynamics

   -  Solvation Energies

Programming Style
~~~~~~~~~~~~~~~~~

The programming style of DL_POLY_4 is intended to be as uniform as possible. The
following stylistic rules apply throughout. Potential contributors of
code are requested to note the stylistic convention.

Programming Language
--------------------

DL_POLY_4 is written in free format :index:`FORTRAN90`. In DL_POLY_4 we have adopted the convention
of *explicit type declaration* i.e. we have used

::

   Implicit None

in all subroutines. Thus all variables must be given an explicit type:
``Integer, Real( Kind = wp)``, etc.

Modularisation and Intent
-------------------------

DL_POLY_4 exploits the full potential of the modularisation concept in :index:`FORTRAN90`.
Variables having in common description of certain feature or method in DL_POLY_4
are grouped in modules. This simplifies subroutines’ calling sequences
and decreases error-proneness in programming as subroutines must define
what they use and from which module. To decrease error-proneness
further, arguments that are passed in calling sequences of functions or
subroutines have defined intent, i.e. whether they are to be:

-  passed in only (``Intent (In)``) - the argument is not allowed to be
   changed by the routine

-  passed out only (``Intent (Out)``) - the “coming in” value of the
   argument is unimportant

-  passed in both directions in and out (``Intent (InOut)``) - the “coming
   in” value of the argument is important and the argument is allowed to
   be changed.

Memory Management
-----------------

DL_POLY_4 exploits the dynamic array allocation features of :index:`FORTRAN90` to assign
the necessary array dimensions.

Target Platforms
----------------

DL_POLY_4 is intended for distributed memory parallel computers.

Compilation of DL_POLY_4 in parallel mode requires **only** a :index:`FORTRAN90` compiler and
Message Passing Interface (MPI) to handle communications. Compilation of DL_POLY_4
in serial mode is also possible and requires **only** a :index:`FORTRAN90` compiler.

Internal Documentation
----------------------

All subroutines are supplied with a header block of :index:`FORTRAN90` comment
(!) records giving:

#. The name of the author and/or modifying author

#. The version number or date of production

#. A brief description of the function of the subroutine

#. A copyright statement.

Elsewhere :index:`FORTRAN90` comment cards (!) are used liberally.

.. _precision:

FORTRAN90 Parameters and Arithmetic Precision
---------------------------------------------

All global parameters defined by the :index:`FORTRAN90` parameter statements are
specified in the module file: setup_module`, which is included at
compilation time in all subroutines requiring the parameters. All
parameters specified in ``setup_module`` are described by one or more
comment cards.

One super-global parameter is defined at compilation time in the
``kinds_f90`` module file specifying the working precision (``wp``) by kind for
real and complex variables and parameters. The default is 64-bit
(double) precision, i.e. ``Real(wp)``. Users wishing to compile the code
with quadruple precision must ensure that their architecture and
FORTRAN90 compiler can allow that and then change the default in
``kinds_f90``. Changing the precision to anything else that is allowed by
the FORTRAN90 compiler and the machine architecture must also be
compliant with the MPI working precision ``mpi_wp`` as defined in
``comms_module`` (in such cases users must correct for that in there).


.. _units:

Units
-----

Internally all :index:`DL_POLY_4<units;DL_POLY>` subroutines and functions assume the use of the following
defined *molecular units*:

-  The unit of time (:math:`t_{o}`) is :math:`1 \times 10^{-12}` seconds
   (i.e. picoseconds)

-  The unit of length (:math:`\ell_{o}`) is :math:`1 \times 10^{-10}`
   metres (i.e. Ångstroms)

-  The unit of mass (:math:`m_{o}`) is :math:`1.6605402 \times 10^{-27}`
   kilograms (i.e. Daltons - atomic mass units)

-  The unit of charge (:math:`q_{o}`) is
   :math:`1.60217733 \times 10^{-19}` Coulombs (i.e. electrons - units
   of proton charge)

-  The unit of energy (:math:`E_{o}=m_{o}(\ell_{o}/t_{o})^{2}`) is
   :math:`1.6605402 \times 10^{-23}` Joules (10 J mol\ :math:`^{-1}`)

-  The unit of :index:`pressure<units;pressure>` (:math:`{\cal P}_{o}=E_{o}\ell_{o}^{-3}`) is
   :math:`1.6605402 \times 10^{7}` Pascals
   (:math:`163.882576` atmospheres)

-  Planck’s constant (:math:`\hbar`) which is
   :math:`6.350780668 \times E_{o} t_{o}`  .

In addition, the following conversion factors are used:

-  The :index:`coulombic<potential;electrostatics>` conversion factor (:math:`\gamma_{o}`) is:

   .. math:: \gamma_{o} = \frac{1}{E_{o}} \left[ \frac{q_{o}^{2}}{4 \pi \epsilon_{o} \ell_{o}} \right] = 138935.4835~~,

   such that:

   .. math:: U_{\tt MKS}=E_{o}\gamma_{o}U_{\tt Internal}~~,

   where :math:`U` represents the configuration energy.

-  The Boltzmann factor (:math:`k_{B}`) is
   :math:`0.831451115~E_{o}` K:math:`^{-1}`, such that:

   .. math:: T=E_{kin}/k_{B}

   represents the conversion from kinetic energy (in internal units) to
   temperature (in Kelvin).

.. note::
   In the DL_POLY_4 OUTPUT file, the print out of :index:`pressure<units;pressure>` is in units of katms
   (kilo-atmospheres) at all times. The unit of energy is either DL_POLY
   units specified above, or in other units specified by the user at run
   time (see Section :ref:`field-file`. The default is the
   DL_POLY unit.

Externally, DL_POLY_4 accepts information in its own specific formatting as
described in Section :ref:`input-files`. Irrespective of
formatting rules, all values provided to define input entities are read
in DL_POLY units (except otherwise specified as in the case of energy
units) or their composite mixture representing the corresponding entity
physically, i.e. velocities’ components are in Ångstroms/picosecond.

**Exception:** It should be noted that when DL_POLY_4 is used in a DPD mode (see
Section :ref:`dpd` and  :ref:`Appendix A<DPD-all>`) then
the meaning of the molecular units is somewhat lost and it is only the
interrelationship between units that is important (which can be
exploited by the modeller)! The fundamental units for a DPD simulation
are related those of mass :math:`[M]`, length :math:`[L]` and energy
:math:`[E]` - all irrespectively of the actually chosen energy units by
the **UNITS** directive in the FIELD file. Therefore, the DPD unit of time
is equivalent to :math:`[L]\sqrt{[M]/[E]}` while temperature (in the
form :math:`k_{B}T`) is defined as two-thirds of the kinetic energy of
the system’s particles. Similarly, volume is in units of :math:`[L]^{3}`
and pressure in :math:`[E]/[L]^{3}`.

Error Messages
--------------

All errors detected by DL_POLY_4 during run time initiate a call to the subroutine
``error``, which prints an error message in the standard output file and
terminates the program. All terminations of the program are global (i.e.
every node of the parallel computer will be informed of the termination
condition and stop executing).

In addition to terminal error messages, DL_POLY_4 will sometimes print warning
messages. These indicate that the code has detected something that is
unusual or inconsistent. The detection is non-fatal, but the user should
make sure that the warning does represent a harmless condition.


.. _directory-structure:

Directory Structure
~~~~~~~~~~~~~~~~~~~

.. index::
   single: sub-directory;manual
   single: sub-directory;source
   single: sub-directory;build
   single: sub-directory;cmake
   single: sub-directory;utils
   single: sub-directory;execute 
   single: sub-directory;data 
   single: sub-directory;bench 
   single: sub-directory;java
   single: sub-directory;utility

The entire DL_POLY_4 package is stored in a UNIX directory structure. The topmost
directory is named *dl_poly_4.nn*, where *nn* is a generation number.
Beneath this directory are several sub-directories named:
*manual*, *source*, *build*, *cmake*, *utils*, *execute*, 
*data*, *bench*, *java*, and *utility*.

Briefly, the content of each sub-directory is as follows:

.. list-table:: 

   * - sub-directory
     - contents 
   * - *manual* 
     - DL_POLY_4 main user manual and DL_POLY_4 Java GUI manual 
   * - *source* 
     - primary subroutines for the DL_POLY_4 package 
   * - *build* 
     - makefiles to assemble and compile DL_POLY_4 source 
   * - *cmake*
     - contains files needed for DL_POLY_4 ``cmake`` build system
   * - *utils*
     - contains a series of scripts needed for testing 
   * - *execute* 
     - the DL_POLY_4 run-time directory 
   * - *data* 
     - example input and output files for DL_POLY_4 
   * - *bench* 
     -  large test cases suitable for benchmarking 
   * - *java*
     - directory of Java and FORTRAN routines for the Java GUI
   * - *utility*
     - directory of routines donated by DL_POLY_4 users

A more detailed description of each sub-directory follows.

.. index:: single: sub-directory;source

The *source* Sub-directory
--------------------------

In this sub-directory all the essential source code for DL_POLY_4, excluding the
utility software is stored. In keeping with the ‘package’ concept of DL_POLY_4,
it does not contain any complete programs; these are assembled at
compile time using an appropriate makefile. The subroutines in this
sub-directory are documented in Chapter :ref:`source-code_sec`.

.. index:: single: sub-directory;build

The *build* Sub-directory
-------------------------

This sub-directory contains legacy makefiles for the creation (i.e.
compilation and linking) of the DL_POLY_4 simulation program. The makefiles
supplied select the appropriate subroutines from the *source*
sub-directory and deposit the executable program in the *execute*
directory. Building DL_POLY_4 by using these legacy makefiles is described in
Section :ref:`compilation`.

.. index:: single: sub-directory;cmake

The *cmake* Sub-directory
-------------------------

This sub-directory contains necessary scripts and information needed for DL_POLY_4
the CMake system. Building DL_POLY_4 with ``cmake`` is described in
Section :ref:`compilation`.

.. index:: single: sub-directory;utils

The *utils* Sub-directory
-------------------------

This sub-directory contains a framework of scripts needed by DL_POLY_4 developers
for testing purposes. The general user is welcome to look and learn from
it. The scripts are the documentation themselves.

.. index:: single: sub-directory;execute

The *execute* Sub-directory
---------------------------

In the supplied version of DL_POLY_4, this sub-directory contains only a few
macros for copying and storing data from and to the *data* sub-directory
and for submitting programs for execution (see
:ref:`Appendix C<macros>`). However, if the DL_POLY_4 program is assembled
by using a legacy makefile, the executable will be placed in this
sub-directory and could be used from here. Then output files from a job
run in here will also appear here, so users may find it convenient to
use this sub-directory as originally intended. (The experienced user is
not at all required to use DL_POLY_4 this way however.)

.. index:: single: sub-directory;data

The *data* Sub-directory
------------------------

This sub-directory contains examples of input and output files for
testing the released version of DL_POLY_4. The examples of input data are copied
into the *execute* sub-directory when a program is being tested. The test
cases are documented in Chapter :ref:`examples_sec`. Note that these
are no longer within the distribution of any DL_POLY version but are
downloaded when building with testing enabled in cmake.

.. index:: single: sub-directory;java

The *java* Sub-directory
------------------------

The DL_POLY_4 Java Graphical User Interface (:index:`GUI`) is based on the Java\ |reg| language
developed by Sun\ |reg|. The Java\ |reg| source code for this GUI is to be found in
this sub-directory. The source is complete and sufficient to create a
working GUI, provided the user has installed the Java\ |reg| Development Kit,
(1.7 or above) which is available free from Sun\ |reg| at `<http://java.sun.com/>`_.
The GUI, once compiled, may be executed on any machine where Java\ |reg| is
installed :cite:`smith-gui`.

.. index:: single: sub-directory;utility

The *utility* Sub-directory
---------------------------

This sub-directory contains assorted routines donated by DL_POLY users.
Potential users should note that these routines are **unsupported** and come
**without any guarantee or liability whatsoever**. They should be regarded
as potentially useful resources to be hacked into shape as needed by the
user. Some of the various routines in this sub-directory are documented
in the DL_POLY_Classic User Manual. Users who devise their own utilities are advised to
store them in the *utility* sub-directory.


.. _distribution:
.. index:: 
   single: user registration 
   single: DL_POLY_4 software licence

Obtaining the Source Code
~~~~~~~~~~~~~~~~~~~~~~~~~

To obtain a copy of DL_POLY_4 it is necessary to have internet connection. Log on
to the DL_POLY :index:`website<WWW>` - `<http://www.ccp5.ac.uk/DL\_POLY/>`_, and follow the links to the DL_POLY_4 registration
page, where you will firstly be shown the DL_POLY_4 academic software licence (see
:ref:`Appendix E<readme>`), which details the terms and
conditions under which the code will be supplied. **By proceeding further
with the registration and download process you are signalling your
acceptance of the terms of this licence.** Click the ‘Registration’ button
to find the registration page, where you will be invited to enter your
name, address and e-mail address. The code is supplied free of charge to
**academic** users, but **commercial** users will be required to purchase a
software licence.

Once the online registration has been completed, information on
downloading the DL_POLY_4 source code will be sent by e-mail, so **it is therefore
essential to supply a correct e-mail address.**

The *data* and *bench* subdirectories of DL_POLY_4 are not issued in the standard
package, but can be downloaded directly from the FTP site (in the ``ccp5/DL_POLY/DL_POLY_4.0/``
directory).

**Note:** Daresbury Laboratory is the **sole centre** for the distribution of and
copies obtained from elsewhere will be regarded as illegal and will not
be supported.

OS and Hardware Specific Ports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Note that no support is offered for these highly specific developments!**



.. _otherInfoSection:

Other Information
~~~~~~~~~~~~~~~~~

The DL_POLY :index:`website<WWW>` - `<http://www.ccp5.ac.uk/DL\_POLY/>`_, provides additional information in the form of

#. Access to all documentation (including licences)

#. Frequently asked questions

#. Bug reports

#. Access to the DL_Software portal.

Daresbury Laboratory also maintains a associated electronic mailing
list, *dl_poly_4_news*, to which all registered DL_POLY_4 users are automatically
subscribed. It is via this list that error reports and announcements of
new versions are made. If you are a DL_POLY_4 user, but not on this list you may
request to be added by sending a mail message to `majordomo@dl.ac.uk`_ with
the one-line message: :math:`subscribe~dl\_poly\_4\_news`.

The DL_Software **Portal** is a web based centre for all DL_POLY users to
exchange comments and queries. You may access the forum through the
DL_POLY website. A registration (and vetting) process is required before
you can use the forum, but it is open, in principle, to everyone.
