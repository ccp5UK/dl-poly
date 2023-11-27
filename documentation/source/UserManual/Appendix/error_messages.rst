.. _error-messages:

Appendix D: DL_POLY_4 Error Messages & User Action
==================================================

Introduction
------------

.. index:: single: error messages

In this appendix we document the error messages encoded in DL_POLY_4 and
the recommended user action. The correct response is described as the
**standard user response** in the appropriate sections below, to which
the user should refer before acting on the error encountered.

The reader should also be aware that some of the error messages listed
below may be either disabled in, or absent from, the public version of .
Note that the wording of some of the messages may have changed over
time, usually to provide more specific information. The most recent
wording appears below.

The Standard User Response
~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 uses :index:`FORTRAN90` dynamic array allocation to set the array sizes
at run time. This means that a single executable may be compiled to over
all the likely uses of the code. It is not foolproof however. Sometimes
an estimate of the required array sizes is difficult to obtain and the
calculated value may be too small. For this reason DL_POLY_4 retains
array dimension checks and will terminate when an array bound error
occurs.

When a dimension error occurs, the **standard user response** is to
edit the DL_POLY_4
subroutine ``set_bounds``. Locate where the variable defining the
array dimension is fixed and increase accordingly. To do this you
should make use of the dimension information that DL_POLY_4 prints in
the OUTPUT file prior to termination. If no information is supplied,
simply doubling the size of the variable will usually do the trick. If
the variable concerned is defined in one of the support subroutines
scan_config, scan_field, scan_control you will need to insert a new
line in ``set_bounds`` to redefine it - after the relevant subroutine
has been called! Finally the code must be recompiled, as in this case
it will only be necessary to recompile ``set_bounds`` and not the
whole code.

The DL_POLY_4 Error Messages
----------------------------

**Message 1**: error - word_2_real failure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The semantics in some of the INPUT files is wrong. DL_POLY_4 has tried
to read a number but the has found a word in non-number format.

*Action*:

Look into your INPUT files and correct the semantics where appropriate
and resubmit. DL_POLY_4 will have printed out in the OUTPUT file what
the found non-uniform word is.

**Message 2**: error - too many atom types in FIELD (scan_field)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error arises when DL_POLY_4 scans the FIELD file and discovers that
there are too many different types of atoms in the system (i.e. the
number of unique atom types exceeds the 1000).

*Action*:

Increase the number of allowed atom types (mmk) in scan_field, recompile
and resubmit.

**Message 3**: error - unknown directive found in CONTROL file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error most likely arises when a directive is misspelt in the
CONTROL file.

*Action*:

Locate the erroneous directive in the CONTROL file and correct error and
resubmit.

**Message 4**: error - unknown directive found in FIELD file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error most likely arises when a directive is misspelt or is
encountered in an incorrect location in the FIELD file, which can happen
if too few or too many data records are included.

*Action*:

Locate the erroneous directive in the FIELD file and correct error and
resubmit.

**Message 5**: error - unknown energy unit requested
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DL_POLY_4 FIELD file permits a choice of units for input of energy
parameters. These may be: electron-Volts (**eV**); k-calories per mol
(**kcal**/mol); k-Joules per mol (**kJ**/mol); Kelvin per Boltzmann
(**K**\ elvin/Boltzmann); or the DL_POLY_4 internal units, 10 Joules per
mol (**internal**). There is no default value. Failure to specify any of
these correctly, or reference to other energy units, will result in this
error message. See documentation of the FIELD file.

*Action*:

Correct energy keyword on **units** directive in FIELD file and
resubmit.

**Message 6**: error - energy unit not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A **units** directive is mandatory in the FIELD file. This error
indicates that DL_POLY_4 has failed to find the required record.

*Action*:

Add **units** directive to FIELD file and resubmit.

**Message 7**: error - selected external field incompatible with selected ensemble (NVE only!!!)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

Change the external field directive in FIELD file and or the type of
ensemble in CONTROL and resubmit.

**Message 8**: error - ewald precision must be a POSITIVE real number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ewald precision must be a positive non-zero real number. For example
10e-5 is accepted as a standard.

*Action*:

Put a correct number at the "ewald precision" directive in the CONTROL
file and resubmit.

**Message 9**: error - ewald sum parameters must be well defined
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ewald sum parameters must be well defined.

*Action*:

Referer to the manual and references within for understanding the
meaning of the parameters and how to chose them. Alternatively, try
using the “ewald precision” CONTROL directive with a sensible precision
value, of say 10\ :math:`_{-5}`.

**Message 10**: error - too many molecular types specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This indicates an erroneous FIELD file or
corrupted DL_POLY_4 executable. Unlike , DL_POLY_4 does not have a set
limit on the number of kinds of molecules it can handle in any
simulation (this is not the same as the number of molecules).

*Action*:

Examine FIELD for erroneous directives, correct and resubmit.

**Message 11**: error - duplicate molecule directive in FIELD file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The number of different types of molecules in a simulation should only
be specified once. If DL_POLY_4 encounters more than one molecules
directive, it will terminate execution.

*Action*:

Locate the extra **molecule** directive in the FIELD file and remove and
resubmit.

**Message 12**: error - unknown molecule directive in FIELD file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once DL_POLY_4 encounters the **molecules** directive in the FIELD file,
it assumes the following records will supply data describing the
intra-molecular :index:`force field`. It does not then expect to encounter
directives not related to these data. This error message results if it
encounters a unrelated directive. The most probable cause is incomplete
specification of the data (e.g. when the **finish** directive has been
omitted.)

*Action*:

Check the molecular data entries in the FIELD file, correct and
resubmit.

**Message 13**: error - molecule species not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error arises when DL_POLY_4 encounters :index:`non-bonded<potential;non-bonded>` 
force data in the
FIELD file, *before* the molecular species have been specified. Under
these circumstances it cannot assign the data correctly, and therefore
terminates.

*Action*:

Make sure the molecular data appears before the non-bonded forces data
in the FIELD file and resubmit.

**Message 14**: error - too many unique atom types specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This error most likely arises when the FIELD
file or/and DL_POLY_4 executable are corrupted.

*Action*:

Recompile the program and/or recreate the FIELD file afresh. If no
combination of these works, send the problem to us.

**Message 15**: error - duplicate vdw potential specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In processing the FIELD file, DL_POLY_4 keeps a record of the specified
short range pair potentials as they are read in. If it detects that a
given pair potential has been specified before, no attempt at a
resolution of the ambiguity is made and this error message results. See
specification of FIELD file.

*Action*:

Locate the duplication in the FIELD file, rectify and resubmit.

**Message 16**: error - strange exit from FIELD file processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! It simply means that DL_POLY_4 has ceased
processing the FIELD data, but has not reached the end of the file or
encountered a **close** directive. Probable cause: corruption of the
DL_POLY_4 executable or of the FIELD file. We would be interested to
hear of other reasons!

*Action*:

See action notes on message 14 above.

**Message 17**: error - strange exit from CONTROL file processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! It simply means that DL_POLY_4 has ceased
processing the CONTROL data, but has not reached the end of the file or
encountered a **close** directive. Probable cause: corruption of the
DL_POLY_4 executable or of the FIELD file. We would be interested to
hear of other reasons!

*Action*:

Recompile the program and/or recreate the CONTROL file afresh. If no
combination of these works, send the problem to us.

**Message 18**: error - duplicate three-body potential specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has encountered a repeat specification of a :index:`three-body<potential;three-body>`
potential in the FIELD file.

*Action*:

Locate the duplicate entry, remove and resubmit job.

**Message 19**: error - duplicate four-body potential specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 4-body potential has been duplicated in the FIELD file.

*Action*:

Locate the duplicated :index:`four-body<potential;four-body>` potential, remove and resubmit job.

**Message 20**: error - too many molecule sites specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This error most likely arises when the FIELD
file or/and DL_POLY_4 executable are corrupted.

*Action*:

See action notes on message 14 above.

**Message 21**: error - molecule contains more atoms/sites than declared
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The molecule contains more atom/site entries that it declares in the
beginning.

*Action*:

Recreate or correct the erroneous entries in the FIELD file and try
again.

**Message 22**: error - unsuitable radial increment in TABLETABBNDTABANGTABDIHTABINV file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;tabulated

This arises when the tabulated van der Waals potentials presented in the
TABLE file have an increment that is greater than that used to define
the other potentials in the simulation. Ideally, the increment should be
:math:`r_{\rm cut}/(\texttt{mxgrid}-4)`, where :math:`r_{\rm cut}` is
the largest potential cutoff of all supplied ,for the short range
potentials and the domain decomposition link cell size, and ``mxgrid``
is the parameter defining the length of the interpolation arrays. An
increment less than this is permissible however. The same argument holds
for the tabulated intra-molecular interactions that are possibly
supplied via the TABBND, TABANG, TABDIH and TABINV files. All should
have grids sized less than the generic ``mxgrid-4``.

*Action*:

The tables must be recalculated with an appropriate increment.

**Message 23**: error - incompatible FIELD and TABLE file potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error arises when the specification of the short range potentials
is different in the FIELD and TABLE files. This usually means that the
order of specification of the potentials is different. When DL_POLY_4
finds a change in the order of specification, it assumes that the user
has forgotten to enter one.

*Action*:

Check the FIELD and TABLE files. Make sure that you correctly specify
the pair potentials in the FIELD file, indicating which ones are to be
presented in the TABLE file. Then check the TABLE file to make sure all
the :index:`tabulated<potential;tabulated>` potentials are present in 
the order the FIELD file indicates.

**Message 24**: error - end of file encountered in TABLETABBNDTABANGTABDIHTABINV file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This means the TABLETABBNDTABANGTABDIHTABINV file is incomplete in some
way: either by having too few potentials included, or the number of data
points is incorrect.

*Action*:

Examine the TABLE file contents and regenerate it if it appears to be
incomplete. If it look intact, check that the number of data points
specified is what DL_POLY_4 is expecting.

**Message 25**: error - wrong atom type found in CONFIG file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On reading the input file CONFIG, DL_POLY_4 performs a check to ensure
that the atoms specified in the configuration provided are compatible
with the corresponding FIELD file. This message results if they are not
*or the parallel reading wrongly assumed that CONFIG complies with the
DL_POLY_3/4 style*.

*Action*:

The possibility exists that one or both of the CONFIG or FIELD files has
incorrectly specified the atoms in the system. The user must locate the
ambiguity, using the data printed in the OUTPUT file as a guide, and
make the appropriate alteration. If the reason is in the parallel
reading then produce a new CONFIG using a serial reading and continue
working with it.

**Message 26**: error - neutral group option now redundant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 does not have the neutral group option.

*Action*:

Use the Ewald sum option. (It’s better anyway.)

**Message 27**: error - unit’s member indexed outside molecule’s site range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An intra-molecular or intra-molecular alike interaction (topological)
unit has member/site which is given a number outside the scope of the
molecule it is part of.

*Action*:

Find the erroneous entry in FIELD, correct it and try running DL_POLY_4
again.

**Message 28**: error - wrongly indexed atom entries found in CONFIG file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has detected that the atom indices in the CONFIG file do not
form a contnual and/or non-repeating group of indices.

*Action*:

Make sure the CONFIG file is complies with the DL_POLY_4 standards. You
may use the **no index** option in the CONTROL file to override the
crystalographic sites’ reading from the CONFIG file from reading by
index to reading by order of the atom entries with consecutive
incremental indexing. Using this option assumes that the FIELD topology
description matches the crystalographic sites (atoms entries) in the
CONFIG file by order (consecutively).

**Message 30**: error - too many chemical bonds specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This error most likely arises when the FIELD
file or/and DL_POLY_4 executable are corrupted.

*Action*:

See action notes on message 14 above.

**Message 31**: error - too many chemical bonds per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;bond

DL_POLY_4 limits the number of chemical bond units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxbond`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 32**: error - coincidence of particles in core-shell unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a fault in the definition of a core-shell unit in
the FIELD file. The same particle has been assigned to the core and
shell sites.

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 33**: error - coincidence of particles in constraint bond unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a fault in the definition of a constraint bond unit
in the FIELD file. The same particle has been assigned to the both
sites.

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 34**: error - length of constraint bond unit \ :math:`>=` real space cutoff (rcut)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a constraint bond unit length (FIELD) larger than
the real space cutoff (``rcut``) (CONTROL).

*Action*:

Increase cutoff in CONTROL or decrease the constraint bondlength in
FIELD and resubmit. For small system consider using .

**Message 35**: error - coincidence of particles in chemical bond unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a faulty chemical bond in FIELD (defined between the
same particle).

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 36**: error - only one \*bonds\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one bonds entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 38**: error - outgoing transfer buffer size exceeded in metal_ld_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should not usually happen!

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations. Alternatively, increase ``mxbfxp``
parameter in ``set_bounds`` recompile and resubmit. Send the problem to
us if this is persistent.

**Message 39**: error - incoming data transfer size exceeds limit in metal_ld_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 40**: error - too many bond constraints specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

See action notes on message 14 above.

**Message 41**: error - too many bond constraints per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: constraints;bond

DL_POLY_4 limits the number of bond constraint units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxcons`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 42**: error - undefined direction passed to deport_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Send the problem to us.

**Message 43**: error - outgoing transfer buffer size exceeded in deport_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This may happen in extremely non-equilibrium simulations or usually when
the potentials in use do not hold the system stable.

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations. Alternatively, increase ``mxbfdp``
parameter in ``set_bounds`` recompile and resubmit.

**Message 44**: error - incoming data transfer size exceeds limit in deport_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

See action notes on message 43 above.

**Message 45**: error - too many atoms in CONFIG file or per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can happen in circumstances when indeed the CONFIG file has more
atoms listed than defined in FIELD, or when one of the domains (managed
by an MPI process) has higher particle density than the system average
and contains more particles than allowed by the default based on the
system.

*Action*:

Check if CONFIG and FIELD numbers of particles match. Try executing on
various number of processors. Try using the **densvar** option in
CONTROL to increase ``mxatms`` (alternatively, increase it by hand in
``set_bounds`` and recompile) and resubmit. Send the problem to us if
this is persistent.

**Message 46**: error - undefined direction passed to export_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Send the problem to us.

**Message 47**: error - undefined direction passed to metal_ld_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Send the problem to us.

**Message 48**: error - transfer buffer too small in \*_table_read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

Standard user response. Increase ``mxgrid`` parameter in ``set_bounds``
recompile and resubmit.

**Message 49**: error - frozen shell (core-shell) unit specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DL_POLY_4 option to freeze the location of an atom (i.e. hold it
permanently in one position) is not permitted for the shells in
core-shell units.

*Action*:

Remove the frozen atom option from the FIELD file. Consider using a
non-polarisable atom instead.

**Message 50**: error - too many bond angles specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This error most likely arises when the FIELD
file or/and DL_POLY_4 executable are corrupted.

*Action*:

See action notes on message 14 above.

**Message 51**: error - too many bond angles per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;valence angle

DL_POLY_4 limits the number of valence angle units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxangl`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 52**: error - end of FIELD file encountered
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This message results when DL_POLY_4 reaches the end of the FIELD file,
without having read all the data it expects. Probable causes: missing
data or incorrect specification of integers on the various directives.

*Action*:

Check FIELD file for missing or incorrect data, correct and resubmit.

**Message 53**: error - end of CONTROL file encountered
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This message results when DL_POLY_4 reaches the end of the CONTROL file,
without having read all the data it expects. Probable cause: missing
**finish** directive.

*Action*:

Check CONTROL file, correct and resubmit.

**Message 54**: error - outgoing transfer buffer size exceeded in export_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See naction otes on message 38 above.

**Message 55**: error - end of CONFIG file encountered
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error arises when DL_POLY_4 attempts to read more data from the
CONFIG file than is actually present. The probable cause is an incorrect
or absent CONFIG file, but it may be due to the FIELD file being
incompatible in some way with the CONFIG file.

*Action*:

Check contents of CONFIG file. If you are convinced it is correct, check
the FIELD file for inconsistencies.

**Message 56**: error - incoming data transfer size exceeds limit in export_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 57**: error - too many core-shell units specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

See action notes on message 14 above.

**Message 58**: error - number of atoms in system not conserved
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Either and an atom has been lost in transfer between nodes/domains or
your FIELD is ill defined with respect to what is supplied in
CONFIG/HISTORY.

*Action*:

If this error is issued at start before timestep zero in a simulation
then it is either your FIELD file is ill defined or that your CONFIG
file (or the first frame of your HISTRORY being replayed). Check out for
mistyped number or identities of molecules, atoms, etc. in FIELD and for
mangled/blank lines in CONFIG/HISTORY, or a blank line(s) at the end of
CONFIG or missing FOF (End Of File) character in CONFIG. If this error
is issued after timestep zero in a simulation that is not replaying
HISTORY then it is big trouble and you should report that to the
authors. If it is during replaying HISTORY then your HISTORY file has
corrupted frames and you must correct it before trying again.

**Message 59**: error - too many core-shell units per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: polarisation;shell model

DL_POLY_4 limits the number of core-shell units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxshl`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 60**: error - too many dihedral angles specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

See action notes on message 14 above.

**Message 61**: error - too many dihedral angles per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;dihedral

DL_POLY_4 limits the number of dihedral angle units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxdihd`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 62**: error - too many tethered atoms specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

See action notes on message 14 above.

**Message 63**: error - too many tethered atoms per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;tether

DL_POLY_4 limits the number of tethered atoms in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxteth`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 64**: error - incomplete core-shell unit found in build_book_intra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report problem to authors.

**Message 65**: error - too many excluded pairs specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This error arises when DL_POLY_4 is
identifying the atom pairs that cannot have a pair potential between
them, by virtue of being chemically bonded for example (see subroutine
``build_excl_intra``). Some of the working arrays used in this operation
may be exceeded, resulting in termination of the program.

*Action*:

Contact authors.

**Message 66**: error - coincidence of particles in bond angle unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a fault in the definition of a bond angle in the
FIELD file.

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 67**: error - coincidence of particles in dihedral unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a fault in the definition of a dihedral unit in the
FIELD file.

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 68**: error - coincidence of particles in inversion unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a fault in the definition of a inversion unit in the
FIELD file.

*Action*:

Correct the erroneous entry in FIELD and resubmit.

**Message 69**: error - too many link cells required in three_body_forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The number of link cells required for the build up of the Verlet
neighbour list (as in link_cell_pairs) or the calculation of three- &
four-body as well tersoff forces (as in three_body_forces,
four_body_forces, tersoff_body_forces) in the given model exceeds the
number allowed for by the DL_POLY_4 arrays. Probable cause: your system
has expanded unacceptably much to . This may not be physically sensible!

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations.

**Message 70**: error - constraint_quench failure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: constraints;bond

When a simulation with bond constraints is started, DL_POLY_4 attempts
to extract the kinetic energy of the constrained atom-atom bonds arising
from the assignment of initial random velocities. If this procedure
fails, the program will terminate. The likely cause is a badly generated
initial configuration.

*Action*:

Some help may be gained from increasing the cycle limit, by using the
directive **mxshak** in the CONTROL file. You may also consider reducing
the tolerance of the SHAKE iteration using the directive **shake** in
the CONTROL file. However it is probably better to take a good look at
the starting conditions!

**Message 71**: error - too many metal potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 72**: error - too many tersoff potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 73**: error - too many inversion potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 74**: error - unidentified atom in tersoff potential list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered and erroneous entry for
Tersoff potentials in FIELD.

*Action*:

Correct FIELD and resubmit.

**Message 76**: error - duplicate tersoff potential specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered and erroneous entry for
Tersoff potentials in FIELD.

*Action*:

Correct FIELD and resubmit.

**Message 77**: error - too many inversion angles per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;inversion

DL_POLY_4 limits the number of inversion units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxinv`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 79**: error - tersoff potential cutoff undefined
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered and erroneous entry for
Tersoff potentials in FIELD.

*Action*:

Correct FIELD and resubmit.

**Message 80**: error - too many pair potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 81**: error - unidentified atom in pair potential list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered and erroneous entry for vdw or
metal potentials in FIELD or cited TABle file.

*Action*:

Correct FIELD and/or cited TABle file.

**Message 82**: error - calculated pair potential index too large
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! In checking the vdw and metal potentials
specified in the FIELD file DL_POLY_4 calculates a unique integer
indices that henceforth identify every specific potential within the
program. If this index becomes too large, termination of the program
results.

*Action*:

Report to authors.

**Message 83**: error - too many three-body/angles potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 84**: error - unidentified atom in three-body/angles potential list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered and erroneous entry at
three-body or angles definitions in FIELD.

*Action*:

Correct FIELD and resubmit.

**Message 85**: error - required velocities not in CONFIG file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user attempts to start up a DL_POLY_4 simulation with any type of
**restart** directive (see description of CONTROL file,) the program
will expect the CONFIG file to contain atomic velocities as well as
positions. Termination results if these are not present.

*Action*:

Either replace the CONFIG file with one containing the velocities, or if
not available, remove the **restart ...** directive altogether and let
DL_POLY_4 create the velocities for itself.

**Message 86**: error - calculated three-body potential index too large
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;three-body

This should never happen! DL_POLY_4 has a permitted maximum for the
calculated index for any three-body potential in the system (i.e. as
defined in the FIELD file). If there are :math:`m` distinct types of
atom in the system, the index can possibly range from :math:`1` to
:math:`(m^{2}*(m-1))/2`. If the internally calculated index exceeds this
number, this error reports results.

*Action*:

Report to authors.

**Message 88**: error - legend array exceeded in build_book_intra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second dimension of a legend array has been exceeded.

*Action*:

If you have an intra-molecular (like) interaction present in abundance
in your model that you suspect is driving this out of bound error
increase its legend bound value, ``mxf``\ interaction, at the end of
``scan_field``, recompile and resubmit. If the error persists contact
authors.

**Message 89**: error - too many four-body/dihedrals/inversions potentials specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 90**: error - specified tersoff potentials have different types’
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is not allowed! Only one general type of tersoff potential is
allowed in FIELD as there are no mixing rules between different tersoff
potentials!

*Action*:

Correct your model representation in FIELD and try again.

**Message 91**: error - unidentified atom in four-body/dihedrals/inversions potential list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;four-body

The specification of a four-body or dihedrals or inversions potential in
the FIELD file has referenced an atom type that is unknown.

*Action*:

Locate the errant atom type in the four-body/dihedrals/inversions
potential definition in the FIELD file and correct. Make sure this atom
type is specified by an ``atoms`` directive earlier in the file.

**Message 92**: error - specified metal potentials have different types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;metal

The specified metal interactions in the FIELD file are referencing more
than one generic type of metal potentials. Only one such type is allowed
in the system.

*Action*:

Locate the errant metal type in the metal potential definition in the
FIELD file and correct. Make sure only one metal type is specified for
all relevan atom interactions in the file.

**Message 93**: error - PMFs mixing with rigid bodies not allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

Correct FIELD and resubmit.

**Message 95**: error - error - rcut or (rcut+rpad) :math:`>` minimum of all half-cell widths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order for the minimum image convention to work correctly within , it
is necessary to ensure that the major cutoff, plus its possible padding
distance, applied to the pair interactions does not exceed half the
perpendicular width of the simulation cell. (The perpendicular width is
the shortest distance between opposing cell faces.) Termination results
if this is detected. In NVE and NVT simulations this can only happen at
the start of a simulation, but in NPT and N\ :math:`\mat{\sigma}`\ T, it
may occur at any time.

*Action*:

Supply a cutoff that is less than half the cell width. If running
constant pressure calculations, use a cutoff that will accommodate the
fluctuations in the simulation cell. Study the fluctuations in the
OUTPUT file to help you with this.

**Message 96**: error - incorrect atom totals assignments in metal_ld_set_halo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Big trouble. Report to authors.

**Message 97**: error - constraints mixing with rigid bodies not allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

Correct FIELD and resubmit.

**Message 99**: error - cannot have shells as part of a constraint, rigid body or tether
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

Correct FIELD and resubmit.

**Message 100**: error - core-shell unit separation :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This could only happen if FIELD and CONFIG do not match each other or
CONFIG is damaged.

*Action*:

Regenerate CONFIG (and FIELD) and resubmit.

**Message 101**: error - calculated four-body potential index too large
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! DL_POLY_4 has a permitted maximum for the
calculated index for any four-body potential in the system (i.e. as
defined in the FIELD file). If there are :math:`m` distinct types of
atom in the system, the index can possibly range from :math:`1` to
:math:`(m^{2}*(m+1)*(m+2))/6`. If the internally calculated index
exceeds this number, this error report results.

*Action*:

Report to authors.

**Message 102**: error - rcut \ :math:`<` 2*rcter (maximum cutoff for tersoff potentials)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The nature of the Tersoff interaction requires they have at least twice
shorter cutoff than the standard pair interctions (or the major system
cutoff).

*Action*:

Decrease Tersoff cutoffs in FIELD or increase cutoff in CONTROL and
resubmit.

**Message 103**: error - parameter mxlshp exceeded in pass_shared_units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Various algorithms (constraint and core-shell ones) require that
information about ‘shared’ atoms be passed between nodes. If there are
too many such atoms, the arrays holding the information will be exceeded
and DL_POLY_4 will terminate execution.

*Action*:

Use **densvar** option in CONTROL to increase ``mxlshp`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 104**: error - arrays listme and lstout exceeded in pass_shared_units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should not happen! Dimensions of indicated arrays have been
exceeded.

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations.

**Message 105**: error - shake algorithm (constraints_shake) failed to converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: 
    single: algorithm;SHAKE 
    single: constraints;bond 

The SHAKE algorithm for bond constraints is iterative. If the maximum
number of permitted iterations is exceeded, the program terminates.
Possible causes include: a bad starting configuration; too large a time
step used; incorrect :index:`force field` specification; too high a temperature;
inconsistent constraints (over-constraint) etc..

*Action*:

You may try to increase the limit of iteration cycles in the constraint
subroutines by using the directive **mxshak** and/or decrease the
constraint precision by using the directive **shake** in CONTROL. But
the trouble may be much more likely to be cured by careful consideration
of the physical system being simulated. For example, is the system
stressed in some way? Too far from equilibrium?

**Message 106**: error - neighbour list array too small in link_cell_pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Construction of the :index:`Verlet neighbour list` in subroutine
``link_cell_pairs`` non-bonded (pair) force has exceeded the neighbour
list array dimensions.

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations or increase by hand ``mxlist`` in
``set_bounds``.

**Message 107**: error - too many pairs for rdf look up specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! A possible reason is corruption in FIELD
or/and DL_POLY_4 executable.

*Action*:

See action notes on message 14 above.

**Message 108**: error - unidentified atom in rdf look up list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During reading of RDF look up pairs in FIELD DL_POLY_4 has found an
unlisted previously atom type.

*Action*:

Correct FIELD by either defining the new atom type or changing it to an
already defined one in the erroneous line. Resubmit.

**Message 109**: error - calculated pair rdf index too large
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! In checking the RDF pairs specified in the
FIELD file DL_POLY_4 calculates a unique integer index that henceforth
identify every RDF pair within the program. If this index becomes too
large, termination of the program results.

*Action*:

Report to authors.

**Message 108**: error - duplicate rdf look up pair specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During reading of RDF look up pairs in FIELD DL_POLY_4 has found a
duplicate entry in the list.

*Action*:

Delete the duplicate line and resubmit.

**Message 111**: error - bond constraint unit separation :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! DL_POLY_4 has not been able to find an atom in
a processor domain or its bordering neighbours.

*Action*:

Probable cause: link cells too small. Use larger potential cutoff.
Contact DL_POLY_4 authors.

**Message 112**: error - only one \*constraints\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one constraints entry per molecule in
FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 113**: error - intra-molecular bookkeeping arrays exceeded in deport_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One or more bookkeeping arrays for site-related interactions have been
exceeded.

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations. Alternatively, you will need to print extra
diagnostic data from the ``deport_atomic_data`` subroutine to find which
boded-like contribution has exceeded its assumed limit and then correct
for it in ``set_bounds``, recompile and resubmit.

**Message 114**: error - legend array exceeded in deport_atomic_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The array ``legend`` has been exceeded.

*Action*:

Try increasing parameter ``mxfix`` in ``set_bounds``, recompile and
resubmit. Contact DL_POLY_4 authors if the problem persists.

**Message 115**: error - transfer buffer exceeded in update_shared_units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transfer buffer has been exceeded.

*Action*:

Consider increasing parameter ``mxbfsh`` in ``set_bounds``, recompile
and resubmit. Contact DL_POLY_4 authors if the problem persists.

**Message 116**: error - incorrect atom transfer in update_shared_units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An atom has become misplaced during transfer between nodes.

*Action*:

This happens when the simulation is very numerically unstable. Consider
carefully the physical grounds of your simulation, i.e. are you using
the adiabatic shell model for accounting polarisation with too big a
timestep or too large control distances for the variable timestep, is
the ensemble type NPT or N\ :math:`\mat{\sigma}`\ T and the system
target temperature too close to the melting temperature?

**Message 118**: error - construction error in pass_shared_units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should not happen.

*Action*:

Report to authors.

**Message 120**: error - invalid determinant in matrix inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 occasionally needs to calculate matrix inverses (usually the
inverse of the matrix of cell vectors, which is of size 3 :math:`\times`
3). For safety’s sake a check on the determinant is made, to prevent
inadvertent use of a singular matrix.

*Action*:

Locate the incorrect matrix and fix it - e.g. are cell vectors correct?

**Message 122**: error - FIELD file not found
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 failed to find a FIELD file in your directory.

*Action*:

Supply a valid FIELD file before you start a simulation

**Message 124**: error - CONFIG file not found
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 failed to find a CONFIG file in your directory.

*Action*:

Supply a valid CONFIG file before you start a simulation

**Message 126**: error - CONTROL file not found
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 failed to find a CONTROL file in your directory.

*Action*:

Supply a valid CONTROL file before you start a simulation

**Message 128**: error - chemical bond unit separation :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This could only happen if FIELD and CONFIG do not match each other or if
the instantaneous configuration is ill defined because of generation of
large forces on bonded particles. This may be due to having a badly
defined force-field and/or starting form a configuration which is too
much away from equilibrium.

*Action*:

Regenerate CONFIG (and FIELD) and resubmit. Try topology verification by
using ``nfold 1 1 1`` in CONTROL. Try using options as ``scale``,
``cap``, ``zero`` and ``optimise``. Try using smaller SHAKE tolerance if
constraints are present in the system. You may as well try using the
``variable timestep`` option.

**Message 130**: error - bond angle unit diameter :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See action notes on message 128 above.

*Action*:

See action notes on message 128 above.

**Message 132**: error - dihedral angle unit diameter :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 128 above.

*Action*:

See action notes on message 128 above.

**Message 134**: error - inversion angle unit diameter :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 128 above.

*Action*:

See action notes on message 128 above.

**Message 138**: error - incorrect atom totals assignments in refresh_halo_positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen although, sometimes, it could due to ill
defined force field and/or and/or starting form a configuration which is
too much away from equilibrium.

*Action*:

Try using the ``variable timestep`` option and/or running in serial to
determine if particles gain too much speed and leave domains.

**Message 141**: error - duplicate metal potential specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During reading of metal potentials (pairs of atom types) in FIELD
DL_POLY_4 has found a duplicate pair of atoms in the list.

*Action*:

Delete one of the duplicate entries and resubmit.

**Message 150**: error - unknown van der waals potential selected
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 checks when constructing the interpolation tables for the
short ranged potentials that the potential function requested is one
which is of a form known to the program. If the requested potential form
is unknown, termination of the program results. The most probable cause
of this is the incorrect choice of the potential keyword in the FIELD
file.

*Action*:

Read the DL_POLY_4 documentation and find the potential keyword for the
potential desired.

**Message 151**: error - unknown EAM keyword in TABEAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 checks when constructing the interpolation tables for the EAM
metal potentials that the potential function requested is one which is
of a form known to the program. If the requested potential form is
unknown, termination of the program results. The most probable cause of
this is the incorrect choice of the potential keyword in the FIELD file.

*Action*:

Read the DL_POLY_4 documentation and find the potential keyword for the
potential desired.

**Message 152**: error - undefined direction passed to dpd_v_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Report to authors.

**Message 154**: error - outgoing transfer buffer size exceeded in dpd_v_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 156**: error - incoming data transfer size exceeds limit in dpd_v_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 158**: error - incorrect atom totals assignments in dpd_v_set_halo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Big trouble. Report to authors.

**Message 160**: error - undefined direction passed to statistics_connect_spread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

**Message 163**: error - outgoing transfer buffer size exceeded in statistics_connect_spread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transfer buffer has been exceeded.

*Action*:

Consider using ``densvar`` option in CONTROL for extremely
non-equilibrium simulations. Alternatively, increase ``mxbfss``
parameters in ``set_bounds`` recompile and resubmit.

**Message 164**: error - incoming data transfer size exceeds limit in statistics_connect_spread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 163 above.

*Action*:

See action notes on message 163 above.

**Message 170**: error - too many variables for statistics array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error means the statistics arrays appearing in subroutine
``statistics_collect`` are too small. This should never happen!

*Action*:

Contact DL_POLY_4 authors.

**Message 172**: error - duplicate intra-molecular entries specified
in TABBNDTABANGTABDIHTABINV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A duplicate entry has been encountered in the intra-molecular table
file.

*Action*:

Contact DL_POLY_4 authors.

**Message 200**: error - rdf/z-density buffer array too small in system_revive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error indicates that a global summation buffer array in subroutine
``system_revive`` is too small, i.e mxbuff :math:`<` ``mxgrdf``. This
should never happen!

*Action*:

Contact DL_POLY_4 authors.

**Message 210**: error - only one \*angles\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one angles entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 220**: error - only one \*dihedrals\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one dihedrals entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 230**: error - only one \*inversions\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one inversions entry per molecule in
FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 240**: error - only one \*tethers\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one tethers entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 300**: error - incorrect boundary condition for link-cell algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of link cells in DL_POLY_4 implies the use of appropriate
boundary conditions. This error results if the user specifies octahedral
or dodecahedral boundary conditions, which are only available in .

*Action*:

Correct your boundary condition or consider using .

**Message 305**: error - too few link cells per dimension for many-body and tersoff forces subroutines.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The link cells algorithms for many-body and tersoff forces in cannot
work with less than 3 (secondary) link cells per dimension. This depends
on the cell size widths (as supplied in CONFIG) and the largest system
cut-off (as specified in CONTROL although it may be drawn or overridden
by cutoffs specified as part of some potentials’ parameter sets in
FIELD).

*Action*:

Decrease many-body and tersoff potentials cutoffs or/and number of nodes
or/and increase system size.

**Message 307**: error - link cell algorithm violation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 does not like what you are asking it to do. Probable cause:
the cutoff is too large to use link cells in this case.

*Action*:

Rethink the simulation model; reduce the cutoff or/and number of nodes
or/and increase system size.

**Message 308**: error - link cell algorithm in contention with SPME sum precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 does not like what you are asking it to do. Probable cause:
you ask for SPME precision that is not achievable by the current
settings of the link cell algorithm.

*Action*:

Rethink the simulation model; reduce number of nodes or/and SPME sum
precision or/and increase cutoff.

**Message 340**: error - invalid integration option requested
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has detected an incompatibility in the simulation
instructions, namely that the requested integration algorithm is not
compatible with the physical model. It *may* be possible to override
this error trap, but it is up to the user to establish if this is
sensible.

*Action*:

This is a non-recoverable error, unless the user chooses to override the
restriction.

**Message 350**: error - too few degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error can arise if a small system is being simulated and the number
of constraints applied is too large.

*Action*:

Simulate a larger system or reduce the number of constraints.

**Message 360**: error - degrees of freedom distribution problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen for a dynamically sensical system. This error
arises if a model system contains one or more free, zero mass particles.
Zero mass (mass-less) particles/sites are only allowed for shells in
core-shell units and as part of rigid bodies (mass-less but charged RB
sites).

*Action*:

Inspect your FIELD to find and correct the erroneous entries, and try
again.

**Message 380**: error - simulation temperature not specified or :math:`< 1` K
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to find a **temp** directive in the CONTROL file.

*Action*:

Place a **temp** directive in the CONTROL file, with the required
temperature specified.

**Message 381**: error - simulation timestep not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to find a **timestep** directive in the CONTROL
file.

*Action*:

Place a **timestep** directive in the CONTROL file, with the required
timestep specified.

**Message 382**: error - simulation cutoff not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to find a **cutoff** directive in the CONTROL file.

*Action*:

Place a **cutoff** directive in the CONTROL file, with the required
forces cutoff specified.

**Message 387**: error - system pressure not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The target system pressure has not been specified in the CONTROL file.
Applies to NPT simulations only.

*Action*:

Insert a **press** directive in the CONTROL file specifying the required
system pressure.

**Message 390**: error - npt/nst ensemble requested in non-periodic system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A non-periodic system has no defined volume, hence the NPT algorithm
cannot be applied.

*Action*:

Either simulate the system with a periodic boundary, or use another
ensemble.

**Message 402**: error - van der waals not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;van der Waals

The user has not set any cutoff in CONTROL, (``rvdw``) - the van der
Waals potentials cutoff is needed in order for DL_POLY_4 to proceed.

*Action*:

Supply a cutoff value for the van der Waals terms in the CONTROL file
using the directive ``rvdw``, and resubmit job.

**Message 410**: error - cell not consistent with image convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation cell vectors appearing in the CONFIG file are not
consistent with the specified image convention.

*Action*:

Locate the variable ``imcon`` in the CONFIG file and correct to suit the
cell vectors.

**Message 414**: error - conflicting ensemble options in CONTROL file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one **ensemble** directive in the CONTROL
file.

*Action*:

Locate extra **ensemble** directives in CONTROL file and remove.

**Message 416**: error - conflicting force options in CONTROL file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found incompatible directives in the CONTROL file
specifying the electrostatic interactions options.

*Action*:

Locate the conflicting directives in the CONTROL file and correct.

**Message 430**: error - integration routine not available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A request for a non-existent :index:`ensemble` has been made or a request with
conflicting options that DL_POLY_4 cannot deal with.

*Action*:

Examine the CONTROL and FIELD files and remove inappropriate
specifications.

**Message 432**: error - undefined tersoff potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows that DL_POLY_4 has encountered an unfamiliar entry for
Tersoff potentials in FIELD.

*Action*:

Correct FIELD and resubmit.

**Message 433**: error - rcut must be specified for the Ewald sum precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: Ewald;summation

When specifying the desired precision for the Ewald sum in the CONTROL
file, it is also necessary to specify the real space cutoff ``rcut``.

*Action*:

Place the **cut** directive *before* the **ewald precision** directive
in the CONTROL file and rerun.

**Message 436**: error - unrecognised ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An unknown ensemble option has been specified in the CONTROL file.

*Action*:

Locate **ensemble** directive in the CONTROL file and amend
appropriately.

**Message 440**: error - undefined angular potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A form of angular potential has been requested which DL_POLY_4 does not
recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is possible. Alternatively, you
may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``angles_forces`` will be
required.

**Message 442**: error - undefined three-body potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A form of three-body potential has been requested which DL_POLY_4 does
not recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``three_body_forces`` will
be required.

**Message 443**: error - undefined four-body potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;four-body

DL_POLY_4 has been requested to process a four-body potential it does
not recognise.

*Action*:

Check the FIELD file and make sure the keyword is correctly defined.
Make sure that subroutine
``three_body_forces`` contains the code necessary to deal with the
requested potential. Add the code required if necessary, by amending
subroutines read_field and ``three_body_forces``.

**Message 444**: error - undefined bond potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;bond

DL_POLY_4 has been requested to process a bond potential it does not
recognise.

*Action*:

Check the FIELD file and make sure the keyword is correctly defined.
Make sure that subroutine ``bonds_forces`` contains the code necessary
to deal with the requested potential. Add the code required if
necessary, by amending subroutines read_field and ``bonds_forces``.

**Message 445**: error - r_14 :math:`>` rcut in dihedrals_forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 1-4 coulombic scaling for a dihedral angle bonding cannot be
performed since the 1-4 distance has exceeded the system short range
interaction cutoff, ``rcut``, in subroutine dihedral_forces.

*Action*:

To prevent this error occurring again increase ``rcut``.

**Message 446**: error - undefined electrostatic key in dihedral_forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;electrostatics

The subroutine ``dihedral_forces`` has been requested to process a form
of electrostatic potential it does not recognise.

*Action*:

The error arises because the integer key ``keyfrc`` has an inappropriate
value (which should not happen in the standard version of ). Check that
the FIELD file correctly specifies the potential. Make sure the version
of ``dihedral_forces`` does contain the potential you are specifying.
Report the error to the authors if these checks are correct.

*Action*:

To prevent this error occurring again increase ``rvdw``.

**Message 447**: error - only one \*shells\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one shells entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 448**: error - undefined dihedral potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;dihedral

A form of dihedral potential has been requested which DL_POLY_4 does not
recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``dihedral_forces`` (and
its variants) will be required.

**Message 449**: error - undefined inversion potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;inversion

A form of inversion potential has been encountered which DL_POLY_4 does
not recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``inversions_forces`` will
be required.

**Message 450**: error - undefined tethering potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;tether

A form of tethering potential has been requested which DL_POLY_4 does
not recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``tethers_forces`` will be
required.

**Message 451**: error - three-body potential cutoff undefined
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: potential;three-body

The cutoff radius for a three-body potential has not been defined in the
FIELD file.

*Action*:

Locate the offending three-body force potential in the FIELD file and
add the required cutoff. Resubmit the job.

**Message 452**: error - undefined vdw potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A form of vdw potential has been requested which DL_POLY_4 does not
recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field``, ``vdw_generate``\ \* and
``dihedrals_14_vdw`` will be required.

**Message 453**: error - four-body potential cutoff undefined
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: potential;four-body

The cutoff radius for a four-body potential has not been defined in the
FIELD file.

*Action*:

Locate the offending four-body force potential in the FIELD file and add
the required cutoff. Resubmit the job.

**Message 454**: error - unknown external field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A form of external field potential has been requested which does not
recognise.

*Action*:

Locate the offending potential in the FIELD file and remove. Replace
with one acceptable to DL_POLY_4 if this is reasonable. Alternatively,
you may consider defining the required potential in the code yourself.
Amendments to subroutines ``read_field`` and ``external_field_apply``
will be required.

**Message 456**: error - external field xpis-ton is applied to a layer with at least one frozen particle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a layer to emulate a piston no particle constituting it must be
frozen.

*Action*:

Locate the offending site(s) in the FIELD file and unfreeze the
particles.

**Message 461**: error - undefined metal potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A form of metal potential has been requested which DL_POLY_4 does not
recognise.

*Action*:

Locate erroneous entry in the FIELD file and correct the potental
interaction to one of the allowed ones for metals in .

**Message 462**: error - thermostat friction constant must be\ :math:`~>~0`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A zero or negative value for the :index:`thermostat` friction constant has been
encountered in the CONTROL file.

*Action*:

Locate the **ensemble** directive in the CONTROL file and assign a
positive value to the time constant.

**Message 463**: error - barostat friction constant must be\ :math:`~>~0`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A zero or negative value for the :index:`barostat` friction constant has been
encountered in the CONTROL file.

*Action*:

Locate the **ensemble** directive in the CONTROL file and assign a
positive value to the time constant.

**Message 464**: error - thermostat relaxation time constant must be\ :math:`~>~0`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A zero or negative value for the :index:`thermostat` relaxation time constant has
been encountered in the CONTROL file.

*Action*:

Locate the **ensemble** directive in the CONTROL file and assign a
positive value to the time constant.

**Message 466**: error - barostat relaxation time constant must be\ :math:`~>~0`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A zero or negative value for the :index:`barostat` relaxation time constant has
been encountered in the CONTROL file.

*Action*:

Locate the **ensemble** directive in the CONTROL file and assign a
positive value to the time constant.

**Message 467**: error - rho must not be zero in valid buckingham potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

User specified vdw type buckingham potential has a non-zero force and
zero rho constants. Only both zero or both non-zero are allowed.

*Action*:

Inspect the FIELD file and change the values in question appropriately.

**Message 468**: error - r0 too large for snm potential with current cutoff
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specified location (r0) of the potential minimum for a shifted n-m
potential exceeds the specified potential cutoff. A potential with the
desired minimum cannot be created.

*Action*:

To obtain a potential with the desired minimum it is necessary to
increase the van der Waals cutoff. Locate the ``rvdw`` directive in the
CONTROL file and reset to a magnitude greater than r0. Alternatively
adjust the value of r0 in the FIELD file. Check that the FIELD file is
correctly formatted.

**Message 470**: error - n \ :math:`<` m in definition of n-m potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specification of a n-m potential in the FIELD file implies that the
exponent m is larger than exponent n. (Not all versions of DL_POLY_4 are
affected by this.)

*Action*:

Locate the n-m potential in the FIELD file and reverse the order of the
exponents. Resubmit the job.

**Message 471**: error - rcut \ :math:`<` 2*rctbp (maximum cutoff for three-body potentials)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cutoff for the pair interactions is smaller than twice that for the
three-body interactions. This is a bookkeeping requirement for .

*Action*:

Either use a smaller three-body cutoff, or a larger pair potential
cutoff.

**Message 472**: error - rcut \ :math:`<`  2*rcfbp (maximum cutoff for four-body potentials)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cutoff for the pair interactions is smaller than twice that for the
four-body interactions. This is a bookkeeping requirement for .

*Action*:

Either use a smaller four-body cutoff, or a larger pair potential
cutoff.

**Message 474**: error - conjugate gradient mimimiser cycle limit exceeded
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conjugate gradient minimiser exceeded the iteration limit (100 for
the relaxed shell model, 1000 for the configuration minimiser).

*Action*:

Decrease the respective convergence criterion. Alternatively, you may
try to increase the limit by hand in ``core_shell_relax`` or in
``minimise_relax`` respectively and recompile. However, it is unlikely
that such measures will cure the problem as it is more likely to lay in
the physical description of the system being simulated. For example, are
the core-shell spring constants well defined? Is the system being too
far from equilibrium?

**Message 476**: error - shells MUST all HAVE either zero or non-zero masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The polarisation of ions is accounted via a core-shell model as the
shell dynamics is either relaxed - shells have no mass, or adiabatic -
all shells have non-zero mass.

*Action*:

Choose which model you would like to use in the simulated system and
adapt the shell masses in FIELD to comply with your choice.

**Message 478**: error - shake algorithms (constraints & pmf) failed to converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your system has both bond and PMF constraints. SHAKE (RATTLE_VV1) is
done by combined application of both bond and PMF constraints SHAKE
(RATTLE_VV1) in an iterative manner until the PMF constraint virial
converges to a constant. No such convergence is achieved.

*Action*:

See action notes on message 515 below.

**Message 480**: error - PMF constraint length :math:`>` minimum of all half-cell widths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specified PMF length has exceeded the minimum of all half-cell
widths.

*Action*:

Specify shorter PMF length or increase MD cell dimensions.

**Message 484**: error - only one potential of mean force permitted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only one potential of mean force is permitted in FIELD.

*Action*:

Correct the erroneous entries in FIELD.

**Message 486**: error - only one of the PMF units is permitted to have frozen atoms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only one of the PMF units is permitted to have frozen atoms.

*Action*:

Correct the erroneous entries in FIELD.

**Message 488**: error - too many PMF constraints per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should not happen.

*Action*:

Is the use of PMF constraints in your system physically sound?

**Message 490**: error - local PMF constraint not found locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should not happen.

*Action*:

Is your system physically sound, is your system equilibrated?

**Message 492**: error - a diameter of a PMF unit :math:`>` minimum of all half cell widths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The diameter of a PMF unit has exceeded the minimum of all half-cell
widths.

*Action*:

Consider the physical concept you are trying to imply in the simulation.
Increase MD cell dimensions.

**Message 494**: error - overconstrained PMF units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PMF units are oveconstrained.

*Action*:

DL_POLY_4 algorithms cannot handle overconstrained PMF units. Decrease
the number of constraints on the PMFs.

**Message 497**: error - pmf_quench failure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

See notes on message 515 below.

**Message 498**: error - shake algorithm (pmf_shake) failed to converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Action*:

See action notes on message 515 below.

**Message 499**: error - rattle algorithm (pmf_rattle) failed to converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 515 below.

*Action*:

See action notes on message 515 below.

**Message 500**: error - PMF unit of zero length is not permitted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PMF unit of zero length is found in FIELD. PMF units are either a single
atom or a group of atoms usually forming a chemical molecule.

*Action*:

Correct the erroneous entries in FIELD.

**Message 501**: error - coincidence of particles in PMF unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A PMF unit must be constituted of non-repeating particles!

*Action*:

Correct the erroneous entries in FIELD.

**Message 502**: error - PMF unit member found to be present more than once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A PMF unit is a group of unique (distingushed) atoms/sites. No
repetition of a site is allowed in a PMF unit.

*Action*:

Correct the erroneous entries in FIELD.

**Message 504**: error - cutoff too large for TABLETABBND file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The requested cutoff exceeds the information in the TABLE file or the
TABBND cutoff is larger than half the system cutoff ``rcut``.

*Action*:

In the case when this is received while reading TABLE, reduce the value
of the vdw cutoff (``rvdw``) in the CONTROL file or reconstruct the
TABLE file. In the case when this is received while reading TABBND then
specify a larger ``rcut`` in CONTROL.

**Message 505**: error - EAM metal densities or pair crossfunctions out of range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resulting densities or pair crossfunctions are not defined in the
TABEAM file.

*Action*:

Recreate a TABEAM file with wider interval of defined densities and pair
cross functions.

**Message 506**: error - EAM or MBPC metal densities out of range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resulting densities are not defined in the TABEAM file if EAM is
used or ill defined due to atoms nearly overlapping when MBPC metal
potential is in use..

*Action*:

Recreate a TABEAM file with wider range of densities.

**Message 507**: error - metal density embedding out of range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the case of EAM type of metal interactions this indicates that the
electron density of a particle in the system has exceeded the limits for
which the embedding function for this particle’s type is defined (as
supplied in TABEAM. In the case of Finnis-Sinclair type of metal
interactions, this indicates that the density has become negative.

*Action*:

Reconsider the physical sanity and validity of the metal interactions in
your system and this type of simulation. You MUST change the
interactions’ parameters and/or the way the physical base of your
investigation is handled in MD terms.

**Message 508**: error - EAM metal interaction entry in TABEAM unspecified in FIELD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specified EAM metal interaction entry found in TABEAM is not
specified in FIELD.

*Action*:

For :math:`N` metal atom types there are :math:`(5N+N^{2})/2` EAM
functions in the TABEAM file. One density (:math:`N`) and one embedding
(:math:`N`) function for each atom type and :math:`(N+N^{2})/2`
cross-interaction functions. Fix the table entries and resubmit.

**Message 509**: error - duplicate entry for a pair interaction detected in TABEAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A duplicate cross-interaction function entry is detected in the TABEAM
file.

*Action*:

Remove all duplicate entries in the TABEAM file and resubmit.

**Message 510**: error - duplicate entry for a density function detected in TABEAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A duplicate density function entry is detected in the TABEAM file.

*Action*:

Remove all duplicate entries in the TABEAM file and resubmit.

**Message 511**: error - duplicate entry for an embedding function detected in TABEAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A duplicate embedding function entry is detected in the TABEAM file.

*Action*:

Remove all duplicate entries in the TABEAM file and resubmit.

**Message 512**: error - non-definable vdw/dpd interactions detected in FIELD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A VDW corss-interaction was uncpecified and recovering it by using a
mixing rule proved impossible due to type difference of the single
species potentials.

*Action*:

Rethink your FIELD file interactions before restarting the job with a
new compatible FIELD and possibly CONTROL file.

**Message 513**: error - particle assigned to non-existent domain in read_config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can only happen if particle coordinates do not match the cell
parameters in CONFIG. Probably, due to negligence or numerical
inaccuracy inaccuracy in generation of big supercell from a small one.

*Action*:

Make sure lattice parameters and particle coordinates marry each other.
Increase accuracy when generating a supercell.

**Message 514**: error - allowed image conventions are: 0, 1, 2, 3 and 6
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found unsupported boundary condition specified in CONFIG.

*Action*:

Correct your boundary condition or consider using .

**Message 515**: error - rattle algorithm (constraints_rattle) failed to converge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: 
    single: algorithm;SHAKE 
    signle: constraints;bond 

The RATTLE algorithm for bond constraints is iterative. If the maximum
number of permitted iterations is exceeded, the program terminates.
Possible causes include: incorrect :index:`force field` specification; too high a
temperature; inconsistent constraints (over-constraint) etc..

*Action*:

You may try to increase the limit of iteration cycles in the constraint
subroutines by using the directive **mxshak** and/or decrease the
constraint precision by using the directive **shake** in CONTROL. But
the trouble may be much more likely to be cured by careful consideration
of the physical system being simulated. For example, is the system
stressed in some way? Too far from equilibrium?

**Message 517**: error - allowed configuration information levels are: 0, 1 and 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found an erroneous configuration information level,
:math:`l: 0 .le. l .le. 2`, (i) for the trajectory option in CONTROL or
(ii) in the header of CONFIG.

*Action*:

Correct the error in CONFIG and rerun.

**Message 518**: error - control distances for variable timestep not intact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found the control distances for the variable timestep
algorithm to be in contention with each other.

*Action*:

``mxdis`` MUST BE\ :math:`~>~2.5 \times` ``mndis``. Correct in CONTROL
and rerun.

**Message 519**: error - REVOLD is incompatible or does not exist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Either REVOLD does not exist or its formatting is incompatible.

*Action*:

Change the ``restart`` option in CONTROL and rerun.

**Message 520**: error - domain decomposition failed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A DL_POLY_4 check during the domain decomposition mapping has been
violated. The number of nodes allowed for imcon = 0 is only 1,2,4 and 8!
The number of nodes allowed for imcon = 6 is restricted to 2 along the z
direction! The number of nodes should not be a prime number since these
are not factorisable/decomposable!

*Action*:

You must ensure DL_POLY_4 execution on a number of processors that
complies with the advise above.

**Message 530**: error - pseudo thermostat thickness MUST comply with: 2 Angs \ :math:`<=` thickness :math:`<` a quarter of the minimum MD cell width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a check violated while reading CONTROL.

*Action*:

Correct accordingly in CONTROL and resubmit.

**Message 540**: error - pseudo thermostat MUST only be used in bulk simulations, i.e. imcon MUST be 1, 2 or 3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found a check violated while reading CONTROL.

*Action*:

Correct accordingly in CONTROL nve or in CONFIG (``imcon``) and
resubmit.

**Message 551**: error - REFERENCE not found !!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``defect`` detection option is used in conjunction with ``restart``
but no REFERENCE file is found.

*Action*:

Supply a REFERENCE configuration.

**Message 552**: error - REFERENCE must contain cell parameters !!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REFERENCE MUST contain cell parameters i.e. image convention MUST be
**imcon** :math:`=` 1, 2, 3 or 6.

*Action*:

Supply a properly formatted REFERENCE configuration.

**Message 553**: error - REFERENCE is inconsistent !!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An atom has been lost in transfer between nodes. This should never
happen!

*Action*:

Big trouble. Report problem to authors immediately.

**Message 554**: error - REFERENCE’s format different from CONFIG’s !!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REFERENCE complies to the same rules as CONFIG with the exception that
image convention MUST be **imcon** :math:`=` 1, 2, 3 or 6.

*Action*:

Supply a properly formatted REFERENCE configuartion.

**Message 555**: error - particle assigned to non-existent domain in defects_read_reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 513 above.

*Action*:

See action notes on message 513 above.

**Message 556**: error - too many atoms in REFERENCE file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 45 above.

*Action*:

See action notes on message 45 above.

**Message 557**: error - undefined direction passed to defects_reference_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 42 above.

*Action*:

See action notes on message 42 above.

**Message 558**: error - outgoing transfer buffer exceeded in defects_reference_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 559**: error - incoming data transfer size exceeds limit in defects_reference_export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 38 above.

*Action*:

See action notes on message 38 above.

**Message 560**: error - rdef found to be :math:`>` half the shortest interatomic distance in REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The defect detection option relies on a cutoff, rdef, to define the
vicinity around a site (defined in REFERENCES) in which a particle can
claim to occupy the site. Evidently, rdef MUST be :math:`<` half the
shortest interatomic distance in REFERENCE.

*Action*:

Decrease the value of rdef at directive **defect** in CONTROL.

**Message 570**: error - unsupported image convention (0) for system expansion option nfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

System expansion is possible only for system with periodicity on their
boundaries.

*Action*:

Change the image convention in CONFIG to any other suitable periodic
boundary condition.

**Message 580**: error - replay (HISTORY) option can only be used for structural property recalculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No structural property has been specified for this option to activate
itself.

*Action*:

In CONTROL specify properties for recalculation (RDFs,z-density
profiles, defect detection) or alternatively remove the option.

**Message 585**: error - end of file encountered in HISTORY file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This means that the HISTORY file is incomplete in some way: Either
should you abort the replay (HISTORY) option or provide a fresh HISTORY
file before restart.

*Action*:

In CONTROL specify properties for recalculation (RDFs,z-density
profiles, defect detection) or alternatively remove the option.

**Message 590**: error - uknown minimisation type, only "force", "energy" and "distance" are recognised
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configuration minimisation can take only these three criteria.

*Action*:

In CONTROL specify the criterion you like followed by the needed
arguments.

**Message 600**: error - "impact" option specified more than once in CONTROL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only one instance of the "impact" option is allowed in CONTROL.

*Action*:

Remove any extra instances of the "impact" option in CONTROL.

**Message 610**: error - "impact" applied on particle that is either frozen, or the shell of a core-shell unit or part of a RB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is the user’s responsibility to ensure that impact is initiated on a
"valid" particle.

*Action*:

In CONTROL remove the "impact" directive or correct the particle
identity in it so that it complies with the requirements.

**Message 620**: error - duplicate or mixed intra-molecular entries specified in FIELD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FIELD parser has detected an inconsistency in the description of
bonding interactions. It is the user’s responsibility to ensure that no
duplicate or mixed-up intra-molecular entries are specified in FIELD.

*Action*:

Look at the preceding warning message in OUTPUT and find out which entry
of what intra-molecular-like interaction is at fault. Correct the
bonding description and try running again.

**Message 625**: error - only one \*rigid\* directive per molecule is allowed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has found more than one rigids entry per molecule in FIELD.

*Action*:

Correct the erroneous part in FIELD and resubmit.

**Message 630**: error - too many rigid body units specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen! This indicates an erroneous FIELD file or
corrupted DL_POLY_4 executable. Unlike , DL_POLY_4 does not have a set
limit on the number of rigid body types it can handle in any simulation
(this is not the same as the total number of RBs in the system or per
domain).

*Action*:

Examine FIELD for erroneous directives, correct and resubmit.

**Message 632**: error - rigid body unit MUST have at least 2 sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is likely to be a corrupted FIELD file.

*Action*:

Examine FIELD for erroneous directives, correct and resubmit.

**Message 634**: error - rigid body unit MUST have at least one non-massless site
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No RB dynamics is possible if all sites of a body are massless as no
rotational inertia can be defined!

*Action*:

Examine FIELD for erroneous directives, correct and resubmit.

**Message 638**: error - coincidence of particles in rigid body unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This indicates a corrupted FIELD file as all members of a RB unit must
be destinguishable from one another.

*Action*:

Examine FIELD for erroneous directives, correct and resubmit.

**Message 640**: error - too many rigid body units per domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 limits the number of :index:`rigid body` units in the system to be
simulated (actually, the number to be processed by each node) and checks
for the violation of this. Termination will result if the condition is
violated.

*Action*:

Use **densvar** option in CONTROL to increase ``mxrgd`` (alternatively,
increase it by hand in ``set_bounds`` and recompile) and resubmit.

**Message 642**: error - rigid body unit diameter :math:`>` rcut (the system cutoff)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 domain decomposition limits the size of a RB to a largest
diagonal :math:`<` system cutoff. I.e. the largest RB type is still
within a linked cell volume.

*Action*:

Increase cutoff.

**Message 644**: error - overconstrained rigid body unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a very unlikely message which usually indicates a corrupted
FIELD file or unphysically overconstrained system.

*Action*:

Decrease constraint on the system. Examine FIELD for erroneous
directives, if any, correct and resubmit.

**Message 646**: error - overconstrained constraint unit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a very unlikely message which usually indicates a corrupted
FIELD file or unphysically overconstrained system.

*Action*:

Decrease constraint on the system. Examine FIELD for erroneous
directives, if any, correct and resubmit.

**Message 648**: error - quaternion setup failed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error indicates that the routine ``q_setup`` has failed in
reproducing all the atomic positions in rigid units from the centre of
mass and quaternion vectors it has calculated.

*Action*:

Check the contents of the CONFIG file. DL_POLY_4 builds its local body
description of a rigid unit type from the *first* occurrence of such a
unit in the CONFIG file. The error most likely occurs because subsequent
occurrences were not sufficiently similar to this reference structure.
If the problem persists increase the value of ``tol`` in ``q_setup`` and
recompile. If problems still persist double the value of ``dettest`` in
``rigid_bodies_setup`` and recompile. If you still encounter problems
contact the authors.

**Message 650**: error - failed to find principal axis system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error indicates that the routine ``rigid_bodies_setup`` has failed
to find the principal axis for a rigid unit.

*Action*:

This is an unlikely error. DL_POLY_4 should correctly handle linear,
planar and 3-dimensional rigid units. There is the remote possibility
that the unit has all of its mass-bearing particles frozen while some of
the massless are not or the unit has just one mass-bearing particle.
Another, more likely, possibility, in case of linear molecules is that
the precision of the coordinates of these linear molecules’
constituentsi, as produced by the user, is not good enough, which leads
DL_POLY_4 to accepting it as non-linear while, in fact, it is and then
failing at the current point. It is quite possible, despite considered
as wrong practice, that the user defined system of linear RBs is, in
fact, generated from a system of CBs (3 per RB) which has not been run
in a high enough SHAKE/RATTLE tolerance accuracy (10-̂8 and higher may be
needed). Check the definition of the rigid unit in the CONFIG file, if
sensible report the error to the authors.

**Message 655**: error - FENE bond breaking failure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A FENE type bond was broken.

*Action*:

Examine FIELD for erroneous directives, if any, correct and resubmit.

**Message 660**: error - bond-length :math:`>` cutoff in TABBND or cutoff for PDF collection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A bond has reached a length beyond the cutoff over which (i) its
interactions are defined in TABBND or (ii) its potential distribution
function is sampled.

*Action*:

If there is a TABBND present, reconstruct TABBND with the original
interaction potentials defined over a larger cutoff and try to run the
system with new TABBND. If bonds PDFs are collected then increase the
PDF bond cutoff value in CONTROL and try to run the system.

**Message 670**: error - insufficient electronic temperature cells for TTM heat diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The number of coarse-grained electronic temperature (CET) cells for the
heat diffusion calculations of the two-temperature model (TTM) in any
direction (x, y or z) is less than the number of coarse-grained ionic
temperature (CIT) cells.

*Action*:

Examine the OUTPUT for the number of ionic temperature cells and modify
the **ttm ncet** directive in the CONTROL file to ensure there are at
least as many electronic temperature cells.

**Message 680**: error - rpad too large for calculation of ionic temperatures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The padding distance applied to pair interactions is the maximum
distance a particle may exist beyond a periodic boundary. In the case of
calculating ionic temperatures for TTM, this distance extends beyond the
extent (in any direction) of the ionic temperature cells and thus this
property cannot be reliably calculated.

*Action*:

Reduce the padding distance for pair interactions with the **rpad**
directive in the CONTROL file.

**Message 681**: error - electronic specific heat not fully specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not enough information is given for electronic specific heat capacity
functions (other than a constant value) used in two-temperature model
calculations.

*Action*:

Ensure that two positive parameters are given with either the **ttm
cetanh** or **ttm celin** directives in the CONTROL file.

**Message 682**: error - thermal conductivity of metal not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No information given for thermal conductivity of metal used in
two-temperature model calculations.

*Action*:

Ensure that a positive parameter is given with either the **ttm
keconst** or **ttm kedrude** directives in the CONTROL file.

**Message 683**: error - thermal diffusivity of non-metal not specified
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No information given for thermal diffusivity of non-metal used in
two-temperature model calculations.

*Action*:

Ensure that a positive parameter is given with the **ttm diff**
directive in the CONTROL file.

**Message 684**: error - cannot find or open thermal conductivity table file (Ke.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No readable file (Ke.dat) is available to read tabulated thermal
conductivities.

*Action*:

Ensure that a readable text file named Ke.dat is available in the
directory where DL_POLY_4 is being run. (This must be supplied if the
**ttm ketab** directive is given in the CONTROL file.)

**Message 685**: error - no data found in thermal conductivity table file (Ke.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No tabulated thermal conductivities can be read from the Ke.dat file.

*Action*:

Ensure that the Ke.dat file is formatted correctly in two columns:
temperature and thermal conductivity. Temperatures should be greater
than or equal to 0 kelvin.

**Message 686**: error - cannot find or open volumetric heat capacity table file (Ce.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No readable file (Ce.dat) is available to read tabulated volumetric heat
capacities.

*Action*:

Ensure that a readable text file named Ce.dat is available in the
directory where DL_POLY_4 is being run. (This must be supplied if the
**ttm cetab** directive is given in the CONTROL file.)

**Message 687**: error - no data found in volumetric heat capacity table file (Ce.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No tabulated volumetric heat capacities can be read from the Ce.dat
file.

*Action*:

Ensure that the Ce.dat file is formatted correctly in two columns:
temperature and volumetric heat capacity (equal to product of specific
heat capacity and density). Temperatures should be greater than or equal
to 0 kelvin.

**Message 688**: error - cannot find or open thermal diffusivity table file (De.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No readable file (Ce.dat) is available to read tabulated volumetric heat
capacities.

*Action*:

Ensure that a readable text file named De.dat is available in the
directory where DL_POLY_4 is being run. (This must be supplied if the
**ttm detab** directive is given in the CONTROL file.)

**Message 689**: error - no data found in thermal diffusivity table file (De.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No tabulated thermal diffusivities can be read from the De.dat file.

*Action*:

Ensure that the De.dat file is formatted correctly in two columns:
temperature and thermal diffusivity. Temperatures should be greater than
or equal to 0 kelvin.

**Message 690**: error - cannot find or open coupling constant table file (g.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No readable file (g.dat) is available to read tabulated electron-phonon
coupling constants.

*Action*:

Ensure that a readable text file named g.dat is available in the
directory where DL_POLY_4 is being run. (This must be supplied if the
**ttm gvar** directive is given in the CONTROL file.)

**Message 691**: error - no data found in coupling constant table file (g.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No tabulated thermal conductivities can be read from the g.dat file.

*Action*:

Ensure that the g.dat file is formatted correctly in two columns:
temperature and electron-phonon coupling constant. Temperatures should
be greater than or equal to 0 kelvin.

**Message 692**: error - end of file encountered in table file (Ke.dat, Ce.dat, De.dat or g.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 encountered the end of one of the tabulated files (Ke.dat,
Ce.dat, De.dat or g.dat) for two-temperature model calculations
prematurely.

*Action*:

Check that the tabulated files are not corrupted or incomplete in some
way.

**Message 693**: error - negative electronic temperature: instability in electronic heat diffusion equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A negative (non-physical) electronic temperature has been obtained
during solution of the thermal diffusion equation used in the
two-temperature model. This is an indication of instability in the
numerical solution of this partial differential equation.

*Action*:

This error should not happen due to careful selection of the timestep
used for the explicit difference solver. In some limited circumstances,
it may be possible to improve the solver stability by decreasing the
value of ``fopttstep`` in ``ttm_thermal_diffusion``, which is used to
scale the timestep.

**Message 694**: error - electronic temperature restart file (DUMP_E) does not exist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DUMP_E file (containing electronic temperatures for restarting
simulations with the two-temperature model) does not exist.

*Action*:

A DUMP_E file should be supplied if the **restart** directive is
included in the CONTROL file.

**Message 695**: error - mismatch in electronic temperature lattice sizes between restart (DUMP_E) and CONTROL files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The electronic temperature lattice size (number of CET cells) given in
the DUMP_E restart file for two-temperature model simulations does not
correspond to the size given in the CONTROL file.

*Action*:

Ensure that the **ttm ncet** directive in the CONTROL file matches up
with the three numbers in the first line of the DUMP_E file.

**Message 696**: error - cannot read electronic temperature restart (DUMP_E) file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The electronic temperature lattice size (number of CET cells) given in
the DUMP_E restart file for two-temperature model simulations does not
correspond to the size given in the CONTROL file.

*Action*:

Check that the DUMP_E file is not corrupted or incomplete in some way.

**Message 1000**: error - working precision mismatch between FORTRAN90 and MPI implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to match the available modes of MPI precision for
real numbers to the defined in sc kinds_f90 :index:`FORTRAN90` working precision
``wp`` for real numbers. ``wp`` is a precompile parameter.

*Action*:

This simply mean that ``wp`` must have been changed from its original
value to something else and the new value is not matched by the
``mpi_wp`` variable in ``comms_module``. It is the user’s responsibility
to ensure that ``wp`` and ``mpi_wp`` are compliant. Make the necessary
corrections to sc kinds_f90 and/or ``comms_module``.

**Message 1001**: error - allocation failure in comms_module :math:`->` gcheck_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to find available memory to allocate an array or
arrays, i.e. there is lack of sufficient memory (per node) on the
execution machine.

*Action*:

This may simply mean that your simulation is too large for the machine
you are running on. Consider this before wasting time trying a fix. Try
using more processing nodes if they are available. If this is not an
option investigate the possibility of increasing the heap size for your
application. Talk to your systems support people for advice on how to do
this.

**Message 1002**: error - deallocation failure in comms_module :math:`->` gcheck_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 has failed to deallocate an array or arrays, i.e. to free
memory that is no longer in use.

*Action*:

Talk to your systems support people for advice on how to manage this.

**Message 1003**: error - allocation failure in comms_module :math:`->` gisum_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1004**: error - deallocation failure in comms_module :math:`->` gisum_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1005**: error - allocation failure in comms_module :math:`->` grsum_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1006**: error - deallocation failure in comms_module :math:`->` grsum_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1007**: error - allocation failure in comms_module :math:`->` gimax_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1008**: error - deallocation failure in comms_module :math:`->` gimax_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1009**: error - allocation failure in comms_module :math:`->` grmax_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1010**: error - deallocation failure in comms_module :math:`->` grmax_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1011**: error - allocation failure in parse_module :math:`->` get_record
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1012**: error - deallocation failure in parse_module :math:`->` get_record
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1013**: error - allocation failure in angles_module :math:`->` allocate_angles_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1014**: error - allocation failure in bonds_module :math:`->` allocate_bonds_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1015**: error - allocation failure in core_shell_module :math:`->`
allocate_core_shell_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1016**: error - allocation failure in statistics_module :math:`->` allocate_statitics_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1017**: error - allocation failure in tethers_module :math:`->` allocate_tethers_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1018**: error - allocation failure in constraints_module :math:`->`
allocate_constraints_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1019**: error - allocation failure in external_field_module :math:`->`
allocate_external_field_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1020**: error - allocation failure in dihedrals_module :math:`->` allocate_dihedrals_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1021**: error - allocation failure in inversions_module :math:`->` allocate_inversion_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1022**: error - allocation failure in vdw_module :math:`->` allocate_vdw_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1023**: error - allocation failure in metal_module :math:`->` allocate_metal_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1024**: error - allocation failure in three_body_module :math:`->`
allocate_three_body_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1025**: error - allocation failure in config_module :math:`->` allocate_config_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1026**: error - allocation failure in site_module :math:`->` allocate_site_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1027**: error - allocation failure in tersoff_module :math:`->` alocate_tersoff_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1028**: error - deallocation failure in angles_module :math:`->` deallocate_angles_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1029**: error - deallocation failure in bonds_module :math:`->` deallocate_bonds_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1030**: error - deallocation failure in core_shell_module :math:`->`
deallocate_core_shell_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1031**: error - deallocation failure in tethers_module :math:`->`
deallocate_tethers_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1032**: error - deallocation failure in constraints_module :math:`->`
deallocate_constraints_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1033**: error - deallocation failure in dihedrals_module :math:`->`
deallocate_dihedrals_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1034**: error - deallocation failure in inversions_module :math:`->`
deallocate_inversions_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1035**: error - allocation failure in defects_module :math:`->` allocate_defects_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1036**: error - allocation failure in pmf_module :math:`->` allocate_pmf_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1037**: error - deallocation failure in pmf_module :math:`->` deallocate_pmf_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1038**: error - allocation failure in minimise_module :math:`->` allocate_minimise_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1039**: error - deallocation failure in minimise_module :math:`->`
deallocate_minimise_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1040**: error - allocation failure in ewald_module :math:`->` ewald_allocate_kall_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1041**: error - allocation failure in langevin_module :math:`->`
langevin_allocate_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1042**: error - allocation failure in rigid_bodies_module :math:`->`
allocate_rigid_bodies_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1043**: error - deallocation failure in rigid_bodies_module :math:`->`
deallocate_rigid_bodies_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1044**: error - allocation failure in comms_module :math:`->` gimin_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1045**: error - deallocation failure in comms_module :math:`->` gimin_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1046**: error - allocation failure in comms_module :math:`->` grmin_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1047**: error - deallocation failure in comms_module :math:`->` grmin_vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1048**: error - error - allocation failure in comms_module :math:`->` grsum_matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1049**: error - deallocation failure in comms_module :math:`->` grsum_matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1050**: error - sorted I/O base communicator not set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Possible corruption if ``io_module``. This should never happen!

*Action*:

Make sure you have a clean copy of , compiled without any suspicious
warning messages. Contact authors if the problem persists.

**Message 1053**: error - sorted I/O allocation error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your I/O buffer (and possibly batch) size is too big.

*Action*:

Decrease the value of the I/O buffer (and possibly batch) size in
CONTROL and restart your job.

**Message 1056**: error - unkown write option given to sorted I/O
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Contact authors if the problem persists.

**Message 1059**: error - unknown write level given to sorted I/O
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This should never happen!

*Action*:

Contact authors if the problem persists.

**Message 1060**: error - allocation failure in statistics_module :math:`->` allocate_statitics_connect
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1061**: error - allocation failure in statistics_module :math:`->` deallocate_statitics_connect
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1063**: error - allocation failure in vdw_module :math:`->` allocate_vdw_table_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1066**: error - allocation failure in vdw_module :math:`->` allocate_vdw_direct_fs_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1069**: error - allocation failure in metal_module :math:`->` allocate_metal_table_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1070**: error - allocation failure in ewald_module :math:`->` ewald_allocate_kfrz_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1072**: error - allocation failure in bonds_module :math:`->` allocate_bond_pot_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1073**: error - allocation failure in bonds_module :math:`->` allocate_bond_dst_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1074**: error - allocation failure in angles_module :math:`->` allocate_angl_pot_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1075**: error - allocation failure in angles_module :math:`->` allocate_angl_dst_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1076**: error - allocation failure in dihedrals_module :math:`->` allocate_dihd_pot_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1077**: error - allocation failure in dihedrals_module :math:`->` allocate_dihd_dst_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1078**: error - allocation failure in inversions_module :math:`->` allocate_invr_pot_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1079**: error - allocation failure in inversions_module :math:`->` allocate_invr_dst_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1080**: error - allocation failure in greenkubo_module :math:`->` allocate_greenkubo_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1081**: error - allocation failure in dpd_module :math:`->` allocate_dpd_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1083**: error - allocation failure in ttm_module :math:`->` allocate_ttm_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1084**: error - deallocation failure in ttm_module :math:`->` deallocate_ttm_arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1085**: error - allocation failure in ttm_ion_temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1086**: error - deallocation failure in ttm_ion_temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1087**: error - allocation failure in ttm_thermal_diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1088**: error - deallocation failure in ttm_thermal_diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.

**Message 1089**: error - allocation failure in ttm_track_module :math:`->` depoinit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1001 above.

*Action*:

See action notes on message 1001 above.

**Message 1090**: error - deallocation failure in ttm_track_module :math:`->` depoevolve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See notes on message 1002 above.

*Action*:

See action notes on message 1002 above.
