A Guide to Preparing Input Files
================================

The CONFIG file and the FIELD file can be quite large and unwieldy
particularly if a polymer or biological molecule is involved in the
simulation. This section outlines the paths to follow when trying to
construct files for such systems. The description of the force field in
Chapter :ref:`Force Fields<force-field>` is essential reading. The
various utility routines mentioned in this section are described in
greater detail in the User Manual. Many of these have been incorporated
into the DL_POLY :index:`GUI` :cite:`smith-gui` and may be
conveniently used from there.

Inorganic Materials
-------------------

The utility ``genlat`` can be used to construct the CONFIG file for
relatively simple lattice structures. Input is interactive. The FIELD
file for such systems are normally small and can be constructed by hand.
Otherwise, the input of force field data for crystalline systems is
particularly simple, if no angular forces are required (notable
exceptions to this are zeolites and silicate glasses - see below). Such
systems require only the specification of the atomic types and the
necessary pair forces. The reader is referred to the description of the
DL_POLY_4 FIELD file for further details
(Section:ref:`The FIELD File<field-file>`).

DL_POLY_4 can simulate zeolites and silicate (or other) glasses. Both
these materials require the use of angular forces to describe the local
structure correctly. In both cases the angular terms are included as
*three-body terms*, the forms of which are described in
Chapter :ref:`Force Fields<force-field>`. These terms are entered into
the FIELD file with the pair potentials.

An alternative way of handling zeolites is to treat the zeolite
framework as a kind of macromolecule (see below). Specifying all this is
tedious and is best done computationally: what is required is to
determine the nearest image neighbours of all atoms and assign
appropriate bond and valence angle potentials. What must be avoided at
all costs is specifying the angle potentials *without* specifying :index:`bond<potential;bond>`
potentials. In this case will automatically cancel the :index:`non-bonded<potential;non-bonded>` forces
between atoms linked via :index:`valence<potential;valence angle>` angles and the system will collapse.
The advantage of this method is that the calculation is likely to be
faster than using :index:`three-body<potential;three-body>` forces. This method is not recommended for
amorphous systems.

Macromolecules
--------------

.. index:: WWW

To set up force fields for macromolecules, or indeed any covalent
molecules, it is best to use DL_FIELD - `<http://www.ccp5.ac.uk/DL\_FIELD/>`_. It is a program application tool
developed to facilitate the construction of force field models for
biological molecules and other molecules with complex geometries. For
instance proteins, carbohydrates, polymers and networked molecules such
as graphenes and organic cages. **Although created to assist , is a
separate program suite that requires separate registration!**

The primary functions of DL-FIELD are as follows:

#. **Force field model converter:** DL_FIELD converts the user’s atom models,
   supplied in PDB file format, into input files that are recognisable
   and ready to run with and DL_POLY_4 programs with minimum user’s
   intervention. This basically involves the conversion of the user’s
   atomic configuration in simple xyz coordinates into identifiable atom
   types base on a particular user-selectable potential schemes and then
   automatically generate the DL_POLY configuration file (CONFIG), the
   force field file (FIELD) and a generic control file (CONTROL).

#. **Force field editor:** DL_FIELD allows the user to edit or modify parameters
   of a particular force field scheme in order to produce a customised
   scheme that is specific to a particular simulation model. In
   addition, the standard force field model framework can also be easily
   modified. For instance, introduction of pseudo points and rigid body
   implementation to an otherwise standard potential scheme such as
   CHARMM or AMBER, etc.

#. **Force field library repertoire:** DL_FIELD contains a range of popular
   potential schemes (see below), all described in a single DL_FIELD format that
   are also easily recognisable by the user for maintenance purposes.
   Users can easily expand the existing library to include other new
   molecules.

Force Field Schemes
~~~~~~~~~~~~~~~~~~~

The available force field schemes are as follows:

CHARMM - proteins, ethers, some lipids and carbohydrates.

AMBER - proteins and Glycam for carbohydrates.

OPLSAA - proteins

DREIDING - General force field for covalent molecules.

PCFF - Polyorganics and other covalent molecules.

Model Construction
~~~~~~~~~~~~~~~~~~

DL_FIELD does not have feature to construct molecular models. This can be
achieved by either using DL_POLY :index:`GUI` :cite:`smith-gui` or
any other standard molecular building packages. The output files must be
converted into the PDB format. In the case of proteins, these structures
are usually obtained from data banks such as PDB. These ‘*raw* PBD’
files must first be preprocessed by the user before they are in a
‘*readable* format’ for DL_FIELD. To ensure this, it is advisable that users
take into consideration the following steps:

#. Decide on inclusion/exclusion of and if necessary manually delete
   molecular residues that involve multiple occupancies in crystalline
   structures.

#. Usually, hydrogen atoms are not assigned in the raw PDB file. The
   molecules must therefore be pre-filled with hydrogen atoms
   (protonated) by using any standard packages available. The user must
   ensure that proper care is taken of terminal residues which must also
   be appropriately terminated.

#. Decide on the various charge states of some amino acids, such as
   histidine (HIS), lysine (LYS), glutamic acid (GLU), etc., by adding
   or deleting the appropriate hydrogen atoms. Force field schemes such
   as CHARMM will have different three-letter notations for amino acids
   of different charge states, DL_FIELD will automatically identify these
   differences and assign the appropriate potential parameters
   accordingly.

#. For cysteine (CYS) molecules with disulphide bonds, thiolate hydrogen
   atoms must be removed. DL_FIELD will automatically define disulphide bonds
   between the molecules, provided the S-S distance is within a
   ‘sensible’ value.

#. DL_FIELD does not solvate the molecules and it is the user’s responsibility to
   add water by using any standard package available (for example the
   DL_POLY GUI :cite:`smith-gui`).

Fore more details or further information, please consult the DL_FIELD manual
and :index:`website<WWW>`. `<http://www.ccp5.ac.uk/DL\_FIELD/>`_

Adding Solvent to a Structure
-----------------------------

The utility ``wateradd`` adds water from an equilibrated configuration
of 256 SPC water molecules at 300 K to fill out the MD cell. The utility
``solvadd`` fills out the MD box with single-site solvent molecules from
a fcc lattice. The FIELD files will then need to be edited to account
for the solvent molecules added to the file.

Hint: to save yourself some work in entering the non-bonded interactions
variables involving solvent sites to the FIELD file put two bogus atoms
of each solvent type at the end of the CONNECT_DAT file (for :index:`AMBER`
:index:`force field`s) the utility ``ambforce`` will then evaluate all the
:index:`non-bonded<potential;non-bonded>` variables required by . Remember to delete the bogus entries
from the CONFIG file before running .

Analysing Results
-----------------

DL_POLY_4 is not designed to calculate every conceivable property you
might wish from a simulation. Apart from some obvious thermodynamic
quantities and radial distribution functions, it does not calculate
anything beyond the atomic trajectories. You must therefore be prepared
to post-process the HISTORY file if you want other information. There
are some utilities in the package to help with this, but the list is far
from exhaustive. In time, we hope to have many more. Our users are
invited to submit code to the DL_POLY_4 *public* library to help with
this.

The utilities available are described in the User Manual. Users should
also be aware that many of these utilities are incorporated into the
DL_POLY GUI :cite:`smith-gui`.
