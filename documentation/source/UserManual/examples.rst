.. _examples_sec:

Examples 
++++++++

Example Simulations
===================

Because of the size of the data files for the DL_POLY_4 standard test
cases, they are not shipped in the standard download of the DL_POLY_4
source. Test files are downloaded automatically when building/ running
CMake with the CMake variable BUILD_TESTING=ON. This can be done as
follows:

::

       folder="build-mpi-testing"
       rm -rf $folder && mkdir $folder && pushd $folder
       cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON 

Unpack the files in the ‘data’ subdirectory using ‘gunzip’ and ‘tar -xf’
to create the ‘TEST_X’ directory.

These are provided to give examples of DL_POLY_4 simulations and
demonstrate a limited set of relevant functionality over a limited
extent of molecular systems’ complexity only. **Without modification,
they are not necessarily appropriate for serious simulation of the given
systems.** In other words, the examples are not warranted to have
well-defined force fields in terms of applicability, transferability and
fullness, nor are they likely to have a well-defined state point (i.e.
initial configurations may be away from equilibrium, if physical at
all).

The README.txt file supplied both in the *data* directory and in the
directory on the CCP5 FTP server provides a list of all example
simulations used as test cases to check that DL_POLY_4 is working
correctly, including those described in more detail below. All the jobs
are of a size suitable to test the code in parallel execution. They may
not be suitable for a single processor computer. The files are stored in
compressed format. The examples can be run by typing
:: 

   select n

from the *execute* directory, where n is the number of the test case.
The select macro will copy the appropriate input files (at least
CONTROL, CONFIG, and FIELD in all cases) to the *execute* directory
ready for execution. The output file OUTPUT may be compared with the
file supplied in the *data* directory.

Example 1: Sodium Chloride
--------------------------

This is a 27,000 ion system with unit electric charges on sodium and
chlorine. Simulation at 500 K with a NVT Berendsen ensemble. The SPME
method is used to calculate the Coulombic interactions.

Example 2: DPMC in Water
------------------------

The system consists of 200 DMPC molecules in 9379 and water molecules.
Simulation at 300 K using NVE ensemble with SPME and RATTLE algorithm
for the constrained motion. The total system size is 51,737 atoms.

Example 3: KNaSi\ :math:`_{2}`\ O\ :math:`_{5}` - Potassium/Sodium Disilicate Glass
-----------------------------------------------------------------------------------

Potassium Sodium disilicate glass (NaKSi:math:`_{2}`\ O\ :math:`_{5}`)
using two and three-body potentials. Some of the two-body potentials are
read from the TABLE file. Simulation at 1000 K using NVT Nosé-Hoover
ensemble with SPME. Cubic periodic boundaries are in use. The total
system size is 69,120 ions.

Example 4: Gramicidin A Molecules in Water
------------------------------------------

The system consists of 8 gramicidin A molecules in aqueous solution
(32,096 water molecules) with total of 99,120 of atoms. Simulation at
300 K using NPT Berendsen ensemble with SPME and SHAKE/RATTLE algorithm
for the constrained motion.

Example 5: SiC with Tersoff Potentials
--------------------------------------

The system consists of 74,088 atoms. Simulation at 300 K using NPT
Nosé-Hoover ensemble with Tersoff forces and no electrostatics.

Example 6: Cu\ :math:`_{3}`\ Au alloy with Sutton-Chen (metal) Potentials
-------------------------------------------------------------------------

The systems consists of 32,000 atoms. Simulation at 300 K using NVT
Nosé-Hoover ensemble with Sutton-Chen forces and no electrostatics.

Example 7: Lipid Bilayer in Water
---------------------------------

The systems consists of 12,428 atoms. Simulation at 300 K using NVT
Berendsen ensemble with SPME and SHAKE/RATTLE algorithm for the
constrained motion.

Examples 8 and 9: MgO with Adiabatic and with Relaxed Shell Models
------------------------------------------------------------------

These system consist of 8,000 (4,000 shells) charged points. Simulation
at 3000 K using NPT Berendsen ensemble with SPME.

Example 10: Potential of Mean Force on K+ in Water
--------------------------------------------------

The system consists of 13,500 (500 PMFs) atoms. Simulation at 300 K
using NPT Berendsen ensemble with SPME and SHAKE/RATTLE algorithm for
the constrained motion.

Example 11: Cu\ :math:`_{3}`\ Au Alloy with Gupta (metal) Potentials
--------------------------------------------------------------------

The system consists of 32,000 atoms. Simulation at 300 K using NVT
Nosé-Hoover ensemble with Gupta forces and no electrostatics.

Example 12: Cu with EAM (metal) Potential
-----------------------------------------

The system consists of 32,000 atoms. Simulation at 300 K using NPT
Berendsen ensemble with EAM tabulated forces and no electrostatics.

Examples 13 and 14: Al with Analytic and with EAM Tabulated Sutton-Chen (metal) Potentials
------------------------------------------------------------------------------------------

The system consists of 32,000 atoms. Simulation at 300 K using NVT Evans
ensemble with Sutton-Chen forces and no electrostatics.

Examples 15: NiAl Alloy with EAM (metal) Potentials
---------------------------------------------------

The system consists of 27,648 atoms. Simulation at 300 K using NVT Evans
ensemble with EAM tabulated forces and no electrostatics.

Examples 16: Fe with Finnis-Sincair (metal) Potential
-----------------------------------------------------

The system consists of 31,250 atoms. Simulation at 300 K using NPT
Berendsen ensemble with Finnis-Sinclair forces and no electrostatics.

Examples 17: Ni with EAM (metal) Potential
------------------------------------------

The system consists of 32,000 atoms. Simulation at 300 K using NPT
Berendsen ensemble with EAM tabulated forces and no electrostatics.

Examples 18 and 19: SPC IceVII Water with CBs and with RBs
----------------------------------------------------------

The system consists of 11,664 (34,992 atoms) water molecules. Simulation
at 25 K using NVE ensemble with CGM force minimisation and SPME
electrostatics.

Example 20: NaCl Molecules in SPC Water Represented as CBs+RBs
--------------------------------------------------------------

The system consists of 64 NaCl ion pairs with 4,480 water molecules
represented by constraint bonds and 4,416 water molecules represented by
ridig bodies. Totalling 26,816 atoms. Simulation at 295 K using NPT
Berendsen ensemble with CGM energy minimisation and SPME electrostatics.

Example 21: TIP4P Water: RBs with a Massless Charged Site
---------------------------------------------------------

The system consists of 7,263 TIP4P rigid body water molecules totaling
29,052 particles. Simulation at 295 K using NPT Berendsen ensemble with
CGM energy minimisation and SPME electrostatics.

Example 22: Ionic Liquid Dimethylimidazolium Chloride as RBs
------------------------------------------------------------

The system consists of 44,352 ions. Simulation at 400 K using NPT
Berendsen ensemble, using both particle and rigid body dynamics with
SPME electrostatics.

Example 23: Calcite Nano-Particles in TIP3P Water
-------------------------------------------------

In this case 600 molecules of calcium carbonate in the calcite structure
form 8 nano-particles which are suspended in 6,904 water molecules,
represented by a flexible 3-centre TIP3P model. Simulation with SPME
electrostatics at 310 K and 1 atmosphere maintained in a Hoover NPT
ensemble. The system consists of 23,712 ions.

Example 24: Iron/Carbon Alloy with 2BEAM (metal) Potentials
-----------------------------------------------------------

In this case a steel alloy of iron and carbon in ratio 35132 to 1651 is
modelled using an EEAM potential forcefield. Simulation at 1000 K and
0 atmosphere is maintained in a Berendsen NPT ensemble. The system
consists of 36,803 particles.

Example 25: Iron/Chromium Alloy with 2BEAM (metal) Potentials
-------------------------------------------------------------

In this case a steel alloy of iron and chromium in ratio 27635 to 4365
is modelled using an 2BEAM potential forcefield. Simulation at 300 K and
0 atmosphere is maintained in an Evans NVT isokinetic ensemble. The
system consists of 32,000 particles.

Examples 26 and 27: Hexane and Methanol Melts, with Full Atomistic and Coarse-Grained Force-Fields
--------------------------------------------------------------------------------------------------

.. index:: single: WWW

These two examples contain a Hexane and a Methanol melt respectively,
(1000 molecules each) modelled by the OPLSAA force-field (FF). Each
system is also supplied in a CG-mapped representation as converted by
VOTCA, `<http://www.votca.org/>`_, or DL_CGMAP
`<http://www.ccp5.ac.uk/projects/ccp5_cg.shtml>`_.

These test cases are to exemplify the Coarse-Graining (CG) procedure
(see Chapter \ :ref:`coarse-graining`), including
FA-to-CG mapping and obtaining the PMF data by means of Boltzmann
Inversion :cite:`reith-03a`. As a result, DL_POLY_4 could be
used for simulating a CG system with numerically defined, tabulated FFs,
see TABBND, TABANG, TABDIH and TABINV files for intra-molecular
potentials, and TABLE for inter-molecular (short-range, VDW) potentials.

Both tests are also available as parts of the tutorial cases from the
VOTCA package :cite:`ruhle-09a`. Therefore, the CONFIG,
CONTROL and FIELD input files are fully consistent with the
corresponding setup files found in the VOTCA tutorial directories
“csg-tutorials/hexane” and “csg-tutorials/methanol’.

Example 28: Butane in CCl\ :math:`_{4}` Solution with Umbrella Sampling via PLUMED
----------------------------------------------------------------------------------

Free Energy calculation for Buthane with respect to the dihedral angle
as collective variable. We use umbrella sampling as implemented in
PLUMED.

PLUMED enabling in CONTROL:

::

   plumed input umbrella.dat

Contents of umbrella.dat:

::

   phi: TORSION ATOMS=1,2,3,4
   restraint-phi: RESTRAINT ARG=phi KAPPA=500 AT=1.20
   PRINT STRIDE=10 ARG=phi,restraint-phi.bias FILE=COLVAR

Two extra output files are generated in this case: OUTPUT.PLUMED and
COLVAR.

.. note::
   
   a DL_POLY_4 version with PLUMED enabled is used for this.


Example 29: Iron with tabulated EAM (metal) Potential, TTM and Cascade
----------------------------------------------------------------------

In this example 54,000 atoms of iron are modelled with a tabulated
embedded-atom potential optimised to produce correct energetics of point
defects and clusters (M07 in :cite:`malerba-10a`). An energy
impact of 10 keV is applied to an atom and the resulting radiation
damage is evolved using the Two-Temperature Model (TTM) to represent
energy transfers due to electron-phonon coupling and electronic stopping
between atoms and a continuum electronic gas
:cite:`zarkadoula-14a`.

This test case produces additional output files: DUMP_E, LATS_E, LATS_I,
PEAK_E and PEAK_I. It also requires an additional input file (Ce.dat) to
supply tabulated heat capacity data required for evolving the electronic
system.

Example 30: Silicon with original Tersoff Potential, TTM and Swift heavy ion irradiation
----------------------------------------------------------------------------------------

This system consists of 200,000 atoms of silicon modelled using an
original Tersoff (T3) potential. The Two-Temperature Model (TTM) is in
use and an energy deposition is applied to the electronic system using a
Gaussian spatial function, an exponentially decaying temporal function
and an electronic stopping power of 50,000 eV/nm. This simulation
represents Swift heavy ion irradiation in silicon, including the
resulting creation of ion tracks :cite:`khara-16a`.

Example 31: Tungsten with extended Finnis-Sinclair Potential, TTM and laser irradiation
---------------------------------------------------------------------------------------

This system consists of 722,672 atoms of tungsten modelled using an
extended Finnis-Sinclair potential. The Two-Temperature Model (TTM) is
in use and an energy deposition is applied to the electronic system
using a spatial function that is homogeneous in x and y directions and
exponentially decaying in the z direction, as well as a Gaussian
temporal function. This energy deposition represents a laser applied to
the surface of a thin film of tungsten :cite:`murphy-15a`
with a surface fluence of 36 mJ/cm\ :math:`^2` and penetration depth of
12.5 nm, causing the film to expand outwards in the z direction.

Additional input files (Ce.dat and g.dat) are required to supply
tabulated heat capacity and electron-phonon coupling values.

Benchmark Cases
===============

DL_POLY_4 benchmark test cases are available to download them from the
CCP5 FTP server as follows:

::

   FTP site : ftp.dl.ac.uk
   Username : anonymous
   Password : your email address
   Directory: ccp5/DL_POLY/DL_POLY_4.0/BENCH

The DL_POLY_4 authors provide these on an "AS IS" terms. For more
information refer to the README.txt file within.
