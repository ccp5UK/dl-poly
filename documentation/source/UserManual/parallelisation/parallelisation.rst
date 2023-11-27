.. _parallelisation:

Parallelisation
===============


DL_POLY_4 is a distributed parallel molecular dynamics package based on
the Domain Decomposition parallelisation strategy
:cite:`todorov-04a,todorov-06a,pinches-91a,rapaport-91b,smith-91a,smith-93a`.
In this section we briefly outline the basic methodology. Users wishing
to add new features DL_POLY_4 will need to be familiar with the
underlying techniques as they are described in the above references.

The Domain Decomposition Strategy
---------------------------------


The Domain Decomposition (DD) strategy
:cite:`todorov-04a,todorov-06a,smith-91a` is one of several
ways to achieve :index:`parallelisation` in MD. Its name derives from the
division of the simulated system into equi-geometrical spatial blocks or
domains, each of which is allocated to a specific processor of a
parallel computer. I.e. the arrays defining the atomic coordinates
:math:`\underline{r}_{i}`, velocities :math:`\underline{v}_{i}` and forces
:math:`\underline{f}_{i}`, for all :math:`N` atoms in the simulated system,
are divided in to sub-arrays of approximate size :math:`N/P`, where
:math:`P` is the number of processors, and allocated to specific
processors. In DL_POLY_4 the domain allocation is handled by the routine
``domains_module`` and the decision of approximate sizes of various
bookkeeping arrays in ``set_bounds``. The division of the configuration
data in this way is based on the location of the atoms in the simulation
cell, such a geometric allocation of system data is the hallmark of DD
algorithms. Note that in order for this strategy to work efficiently,
the simulated system must possess a reasonably uniform density, so that
each processor is allocated almost an equal portion of atom data (as
much as possible). Through this approach the forces computation and
integration of the equations of motion are shared (reasonably) equally
between processors and to a large extent can be computed independently
on each processor. The method is conceptually simple though tricky to
program and is particularly suited to large scale simulations, where
efficiency is highest.

The DD strategy underpinning DL_POLY_4 is based on the link cell
algorithm of Hockney and Eastwood :cite:`hockney-81a` as
implemented by various authors (e.g. Pinches *et al.*
:cite:`pinches-91a` and Rapaport
:cite:`rapaport-91b`). This requires that the cutoff applied
to the interatomic potentials is relatively short ranged. In DL_POLY_4
the link-cell list is build by the routine ``link_cell_pairs``. As with
all DD algorithms, there is a need for the processors to exchange ‘halo
data’, which in the context of link-cells means sending the contents of
the link cells at the boundaries of each domain, to the neighbouring
processors, so that each may have all necessary information to compute
the pair forces acting on the atoms belonging to its allotted domain.
This in DL_POLY_4 is handled by the ``set_halo_particles`` routine.

Systems containing complex molecules present several difficulties. They
often contain ionic species, which usually require Ewald :index:`summation<Ewald;summation>`
methods :cite:`allen-89a,smith-92b`, and *intra*-molecular
interactions in addition to *inter*-molecular forces. Intramolecular
interactions are handled in the same way as in , where each processor is
allocated a subset of intramolecular bonds to deal with. The allocation
in this case is based on the atoms present in the processor’s domain.
The :index:`SHAKE<algorithm;SHAKE>` and :index:`RATTLE<algorithm;RATTLE>` algorithms
:cite:`ryckaert-77a,andersen-83a` require significant
modification. Each processor must deal with the constraint bonds present
in its own domain, but it must also deal with bonds it effectively
shares with its neighbouring processors. This requires each processor to
inform its neighbours whenever it updates the position of a shared atom
during every SHAKE (RATTLE_VV1) cycle (RATTLE_VV2 updates the
velocities), so that all relevant processors may incorporate this update
into its own iterations. In the case of the DD strategy the SHAKE
(RATTLE) algorithm is simpler than for the Replicated Data method of ,
where global updates of the atom positions (merging and splicing) are
required :cite:`smith-93b`. The absence of the merge
requirement means that the DD tailored SHAKE and RATTLE are less
communications dependent and thus more efficient, particularly with
large processor counts.

The DD strategy is applied to complex molecular systems as follows:

#. Using the atomic coordinates :math:`\underline{r}_{i}`, each processor
   calculates the forces acting between the atoms in its domain - this
   requires additional information in the form of the halo data, which
   must be passed from the neighbouring processors beforehand. The
   forces are usually comprised of:

   #. All common forms of :index:`non-bonded<>potential;non-bonded` atom-atom (van der Waals) forces

   #. Atom-atom (and site-site) coulombic :index:`forces<potential;electrostatics>`

   #. Metal-metal (local density dependent) :index:`forces<potential;metal>`

   #. :index:`Tersoff<potential;Tersoff>` (local density dependent) forces (for hydro-carbons)
      :cite:`tersoff-89a`

   #. :index:`Three-body<potential;three-body>` valence :index:`angle<potential;valence angle>` and hydrogen :index:`bond<potential;bond>` forces

   #. Four-body inversion :index:`forces<potential;four-body>`

   #. Ion core-shell :index:`polarisation<polarisation;shell model>`

   #. Tether :index:`forces<potential;tether>`

   #. Chemical bond :index:`forces<potential;chemical bond>`

   #. Valence angle :index:`forces<potential;valence angle>`

   #. Dihedral angle (and improper dihedral angle) :index:`forces<potential;dihedral>`

   #. Inversion angle :index:`forces<potential;inversion>`

   #. External field :index:`forces<potential;external field>`.

#. The computed forces are accumulated in atomic force arrays
   :math:`\underline{f}_{i}` independently on each processor

#. The force arrays are used to update the atomic velocities and
   positions of all the atoms in the domain

#. Any atom which effectively moves from one domain to another, is
   relocated to the neighbouring processor responsible for that domain.

It is important to note that load balancing (i.e. equal and concurrent
use of all processors) is an essential requirement of the overall
algorithm. In DL_POLY_4 this is accomplished quite naturally through the
DD partitioning of the simulated system. Note that this will only work
efficiently if the density of the system is reasonably uniform.
``There are no load balancing algorithms in DL_POLY_4 to compensate for a bad density distribution.``

Distributing the Intramolecular Bonded Terms
--------------------------------------------

.. index:: single: parallelisation;intramolecular terms

The intramolecular terms in DL_POLY_4 are managed through bookkeeping
arrays which list all atoms (sites) involved in a particular interaction
and point to the appropriate arrays of parameters that define the
potential. Distribution of the forces calculations is accomplished by
the following scheme:

#. Every atom (site) in the simulated system is assigned a unique
   ‘global’ index number from :math:`1` to :math:`N`.

#. Every processor maintains a list of the local indices of the atoms in
   its domain. (This is the local atom list.)

#. Every processor also maintains a sorted (in ascending order) local
   list of global atom indices of the atoms in its domain. (This is the
   local sorted atom list.)

#. Every :index:`intramolecular<parallelisation;intramolecular terms>` 
   bonded term :math:`U_{type}` in the system has a
   unique index number :math:`i_{type}`: from :math:`1` to
   :math:`N_{type}` where :math:`type` represents a :index:`bond<potential;bond>`, 
   :index:`angle<potential;valence angle>`, :index:`dihedral<potential;diheral>`, 
   or :index:`inversion<potential;inversion>`. Also attached there with unique index numbers
   are :index:`core-shell<polarisation;shell model>` units, bond constraint units, 
   :index:`PMF<constraints;PMF>` constraint units,
   :index:`rigid body` units and :index:`tethered<potential;tether>` atoms, their 
   definition by site rather than by chemical type.

#. On each processor a pointer array
   :math:`key_{type}(n_{type},i_{type})` carries the indices of the
   specific atoms involved in the potential term labelled
   :math:`i_{type}` . The dimension :math:`n_{type}` will be :math:`1`
   if the term represents a :index:`tether<potential;tether>`, :math:`1,~2` for a 
   :index:`core-shell<polarisation;shell model>` unit
   or a bond :index:`constraint<constraints;bond>` unit or a bond, 
   :math:`1,~2,~3` for a :index:`valence<potential;valence angle>`
   angle and :math:`1,~2,~3,~4` for a :index:`dihedral<potential;dihedral>` 
   or an :index:`inversion<potential;inversion>`,
   :math:`1,..,n_{\texttt{PMF~unit}_{1~or~2}}+1` for a PMF 
   :index:`constraint<constraints;PMF>`
   unit, or :math:`-1,~0,~1,..,n_{\texttt{RB~unit}}` for a :index:`rigid body`
   unit.

#. Using the :math:`key` array, each processor can identify the global
   indices of the atoms in the bond term and can use this in conjunction
   with the local sorted atoms list *and a binary search algorithm* to
   find the atoms in local atom list.

#. Using the local atom identity, the potential energy and force can be
   calculated.

It is worth mentioning that although :index:`rigid body` units are not bearing
any potential parameters, their definition requires that their topology
is distributed in the same manner as the rest of the intra-molecular
like interactions.

Note that, at the start of a simulation DL_POLY_4 allocates individual
bonded interactions to specific processors, based on the domains of the
relevant atoms (DL_POLY_4 routine ``build_book_intra``). This means that
each processor does not have to handle every possible bond term to find
those relevant to its domain. Also this allocation is updated as atoms
move from domain to domain i.e. during the *relocation* process that
follows the integration of the equations of motion (DL_POLY_4 routine
relocate_particles). Thus the allocation of bonded terms is effectively
dynamic, changing in response to local changes.

Distributing the Non-bonded Terms
---------------------------------

.. index:: single: potential;non-bonded

DL_POLY_4 calculates the non-bonded pair interactions using the link
cell algorithm due to Hockney and Eastwood
:cite:`hockney-81a`. In this algorithm a relatively short
ranged potential cutoff (:math:`r_{\rm cut}`) is assumed. The simulation
cell is logically divided into so-called link cells, which have a width
not less than (or equal to) the cutoff distance. It is easy to determine
the identities of the atoms in each link cell. When the pair
interactions are calculated it is already known that atom pairs can only
interact if they are in the same link cell, or are in link cells that
share a common face. Thus using the link cell ‘address’ of each atom,
interacting pairs are located easily and efficiently via the ‘link list’
that identifies the atoms in each link cell. So efficient is this
process that the link list can be recreated every time step at
negligible cost.

For reasons, partly historical, the link list is used to construct a
:index:`Verlet<algorithm;Verlet neighbour list>` neighbour list 
:cite:`allen-89a`. The Verlet list
records the indices of all atoms within the cutoff radius
(:math:`r_{\rm cut}`) of a given atom. The use of a neighbour list is
not strictly necessary in the context of link-cells, but it has the
advantage here of allowing a neat solution to the problem of ‘excluded’
pair interactions arising from the intramolecular terms and frozen atoms
(see below).

In , the neighbour :index:`list<algorithm;Verlet neighbour list>` 
is constructed *simultaneously* on each node,
using the DD adaptation of the link cell algorithm to share the total
burden of the work reasonably equally between nodes. Each node is thus
responsible for a unique set of non-bonded interactions and the
neighbour list is therefore different on each node.

A feature in the construction of the :index:`Verlet<algorithm;Verlet>` neighbour list for
macromolecules is the concept of *excluded atoms*, which arises from the
need to exclude certain atom pairs from the overall list. Which atom
pairs need to be excluded is dependent on the precise nature of the
:index:`force field` model, but as a minimum atom pairs linked via extensible
:index:`bonds<potential;chemical bond>` or :index:`constraints<constraints;bond>` 
and atoms (grouped in pairs) linked via :index:`valence<potential;valence angle>`
angles are probable candidates. The assumption behind this requirement
is that atoms that are formally :index:`bonded<potential;bonded>` in a chemical sense, should not
participate in :index:`non-bonded<potential;non-bonded>` 
interactions. (However, this is not a
universal requirement of all :index:`force field`\ s.) The same considerations are
needed in dealing with charged excluded atoms.

The modifications necessary to handle the excluded and frozen atoms are
as follows. A distributed *excluded atoms list* is constructed by the
DL_POLY_4 routine ``build_excl_intra`` at the start of the simulation
and is then used in conjunction with the Verlet neighbour list builder
``link_cell_pairs`` to ensure that excluded interactions are left out of
the pair force calculations. Note that, completely frozen pairs of atoms
are excluded in the same manner. The excluded atoms list is updated
during the atom relocation process described above (DL_POLY_4 routine
``exchange_particles``).

Once the neighbour list has been constructed, each node of the parallel
computer may proceed independently to calculate the pair force
contributions to the atomic forces (see routine two_body_forces).

.. index:: 
   single: potential;non-bonded 
   single: potential;metal 
   single: potential;Tersoff

The potential energy and forces arising from the non-bonded
interactions, as well as metal and Tersoff interactions are calculated
using interpolation tables. These are generated in the following
routines: ``vdw_generate``, ``metal_generate``,
``metal_table_derivatives`` and ``tersoff_generate``.

Modifications for the Ewald Sum
-------------------------------

.. index::
   single: Ewald;summation

For systems with periodic boundary conditions DL_POLY_4 employs the
Ewald Sum to calculate the coulombic interactions (see
Section \ :ref:`SPME`). It should be noted that DL_POLY_4 uses
only the Smoothed Particle Mesh (SPME) form of the Ewald sum.

Calculation of the real space component in DL_POLY_4 employs the
algorithm for the calculation of the non-bonded interactions outlined
above, since the real space interactions are now short ranged
(implemented in ``ewald_real_forces`` routine).

The reciprocal space component is calculated using Fast Fourier
Transform (FFT) scheme of the SMPE method
:cite:`essmann-95a,bush-06a` as discussed in
Section \ :ref:`SPME`. The parallelisation of this scheme is
entirely handled within the DL_POLY_4 by the 3D FFT routine
``parallel_fft``, (using ``gpfa_module``) which is known as the
Daresbury advanced Fourier Transform, due to I.J. Bush
:cite:`bush-00a`. This routine distributes the SPME charge
array over the processors in a manner that is completely commensurate
with the distribution of the configuration data under the DD strategy.
As a consequence the FFT handles all the necessary communication
implicit in a distributed SPME application. The DL_POLY_4 subroutine
``ewald_spme_forces`` perfoms the bulk of the FFT operations and charge
array construction, while ``spme_forces`` calculates the forces.

Other routines required to calculate the Ewald sum include ``ewald_module``, ``ewald_excl_forces``, ``ewald_frzn_forces`` and ``spme_container``.

Metal Potentials
----------------

The simulation of metals (Section :ref:`metal`) by DL_POLY_4
makes use of density dependent potentials. The dependence on the atomic
density presents no difficulty however, as this class of potentials can
be resolved into pair contributions. This permits the use of the
distributed :index:`Verlet neighbour list` as outlined above. DL_POLY_4
implements these potentials in various subroutines with names beginning
with ``metal_``.

Tersoff, Three-Body and Four-Body Potentials
--------------------------------------------

DL_POLY_4 can calculate Tersoff, three-body and four-body interactions.
Although some of these interactions have similar terms to some
intramolecular ones (three-body to the bond angle and four-body to
inversion angle), these are not dealt with in the same way as the normal
:index:`bonded<potential;bonded>` interactions. They are generally very short ranged and are most
effectively calculated using a link-cell scheme
:cite:`hockney-81a`. No reference is made to the :index:`Verlet
neighbour` list nor the excluded atoms list. It follows that atoms
involved these interactions can interact via non-bonded (pair) forces
and ionic forces also. Note that contributions from frozen pairs of
atoms to these potentials are excluded. The calculation of the Tersoff
three-body and four-body terms is distributed over processors on the
basis of the domain of the central atom in them. DL_POLY_4 implements
these potentials in the following routines ``tersoff_forces``,
``tersoff_generate``, ``three_body_forces`` and ``four_body_forces``.

Globally Summed Properties
--------------------------

The final stage in the DD strategy, is the global summation of different
(by terms of potentials) contributions to energy, virial and stress,
which must be obtained as a global sum of the contributing terms
calculated on all nodes.

The DD strategy does not require a global summation of the forces,
unlike the Replicated Data method used in , which limits communication
overheads and provides smooth parallelisation to large processor counts.

The Parallel (DD tailored) SHAKE and RATTLE Algorithms
------------------------------------------------------

.. index:: 
   single: algorithm;SHAKE 
   single: algorithm;RATTLE

The essentials of the DD tailored SHAKE and RATTLE algorithms (see
Section \ :ref:`shake-rattle`) are as follows:

#. The bond constraints acting in the simulated system are allocated
   between the processors, based on the location (i.e. domain) of the
   atoms involved.

#. Each processor makes a list of the atoms bonded by 
   :index:`constraints<constraints;bond>` it
   must process. Entries are zero if the atom is not bonded.

#. Each processor passes a copy of the array to the neighbouring
   processors which manage the domains in contact with its own. The
   receiving processor compares the incoming list with its own and keeps
   a record of the shared atoms and the processors which share them.

#. In the first stage of the algorithms, the atoms are updated through
   the usual :index:`Verlet<algorithm;Verlet>` algorithm, without regard to the bond 
   :index:`constraints<constraints;bond>`.

#. In the second (iterative) stage of the algorithms, each processor
   calculates the incremental correction vectors for the 
   :index:`bonded<potential;bond>` atoms in
   its own list of bond :index:`constraints<constraints;bond>`. 
   It then sends specific correction
   vectors to all neighbours that share the same atoms, using the
   information compiled in step 3.

#. When all necessary correction vectors have been received and added
   the positions of the constrained atoms are corrected.

#. Steps 5 and 6 are repeated until the bond constraints are converged.

#. Finally, the change in the atom positions from the previous time step
   is used to calculate the atomic velocities.

The compilation of the list of constrained atoms on each processor, and
the circulation of the list (items 1 - 3 above) is done at the start of
the simulation, but thereafter it needs only to be done every time a
constraint bond atom is relocated from one processor to another. In this
respect DD-SHAKE and DD-RATTLE resemble every other intramolecular term.

Since the allocation of constraints is based purely on geometric
considerations, it is not practical to arrange for a strict load
balancing of the DD-SHAKE and DD-RATTLE algorithms. For many systems,
however, this deficiency has little practical impact on performance.

The Parallel Rigid Body Implementation
--------------------------------------

The essentials of the DD tailored RB algorithms (see
Section \ :ref:`rigid`) are as follows:

#. Every processor works out a list of all local and halo atoms that are
   qualified as free (zero entry) or as members of a RB (unit entry.

#. The rigid body units in the simulated system are allocated between
   the processors, based on the location (i.e. domain) of the atoms
   involved.

#. Each processor makes a list of the RB and their constituting atoms
   that are fully or partially owned by the processors domain.

#. Each processor passes a copy of the array to the neighbouring
   processors which manage the domains in contact with its own. The
   receiving processor compares the incoming list with its own and keeps
   a record of the shared RBs and RBs’ constituent atoms, and the
   processors which share them. *Note that a RB can be shared between up
   to*\ **eight**\ *domains!*

#. The dynamics of each RB is calculated in full on each domain but
   domains only update :math:`\{\underline{r},\underline{v},\underline{f}\}` of RB atoms
   which they own. *Note that a site/atom belongs to*\ **one and only
   one**\ *domain at a time (no sharing) !*

#. Strict bookkeeping is necessary to avoid multiple counting of kinetic
   properties. :math:`\{\underline{r},\underline{v},\underline{v}\}` updates are necessary
   for halo parts (particles) of partially shared RBs. For all domains
   the kinetic contributions from each fully or partially present RB are
   evaluated in full and then waited with the ratio - number of RB’s
   sites local to the domain to total RB’s sites, and then globally
   summed.

The compilation of the lists in items 1 - 3 above and their circulation
of the list is done at the start of the simulation, but thereafter these
need updating on a local level every time a RB site/atom is relocated
from one processor to another. In this respect RBs topology transfer
resembles every other intramolecular term.

Since the allocation of RBs is based purely on geometric considerations,
it is not practical to arrange for a strict load balancing. For many
systems, however, this deficiency has little practical impact on
performance.
