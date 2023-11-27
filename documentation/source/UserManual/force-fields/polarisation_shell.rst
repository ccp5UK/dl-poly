.. _shell-models:

.. index:: polarisation;shell models

Polarisation Shell Models
=========================

.. index:: single:polarisation

An atom or ion is polarisable if it develops a dipole moment when placed
in an electric field. It is commonly expressed by the equation

.. math:: \frac{1}{4\pi\epsilon_{0}\epsilon} \underline{\mu} = \alpha \underline{E}~~,

where :math:`\underline{\mu}` is the induced dipole and :math:`\underline{E}` is the
electric field. The constant :math:`\alpha` is the polarisability.

.. index:: single: polarisation;shell model

In the *static* shell model, also called core-shell model or known as
the Drude model or Druder oscillator
:cite:`drude-00a,drude-00b`, a polarisable atom is
represented by a massive core (often called a nucleus) and a “massless”
shell (also known as a Druder particle), connected by a harmonic spring,
hereafter called the core-shell unit. The core and shell carry different
electric charges, the sum of which equals the charge on the original
*rigid-ion* atom. There is no electrostatic interaction (i.e.
self-interaction) between the core and shell of the same atom.
Non-coulombic interactions arise from the shell alone.

The core-shell interaction is described by a harmonic spring potential
of the form:

.. math:: U_{spring}(r_{ij})=\frac{1}{2}k_{2} r_{ij}^{2}~~,

However, sometimes an anharmonic spring is used, described by a quartic
form:

.. math:: U_{spring}(r_{ij})=\frac{1}{2}k_{2} r_{ij}^{2}+\frac{1}{4!}k_{4} r_{ij}^{4}.

Normally, in practice, :math:`k_{2}` is much larger than :math:`k_{4}`.

The effect of an external electric field, :math:`\underline{E}` is to separate
the core and shell by a distance

.. math:: \underline{d} = q_{s} \underline{E}/k_{2}~~,

giving rise to a *polarisation* dipole

.. math:: \underline{\mu} = q_{s} \underline{d}~~.

The condition of static equilibrium then gives the polarisability as:

.. math:: \alpha = \frac{1}{4\pi\epsilon_{0}\epsilon} q_{s}^{2}/k_{2}~~, \label{druder}
   :label: druder_eq

where :math:`q_{s}` is the shell charge and :math:`k_{2}` is the force
constant of the harmonic spring.

The calculation of the forces, virial and :index:`stress tensor` in this model is
based on that for a diatomic molecule with charged atoms. The part
coming from the spring potential is similar in spirit as for chemical
bonds, equations :eq:`bondf_eq`-\ :eq:`ast_bonds_eq`, while
the electrostatics is as described in the above section. The
relationship between the kinetic energy and the temperature is different
however, as the core-shell unit is permitted only three translational
degrees of freedom, and the degrees of freedom corresponding to rotation
and vibration of the unit are discounted as if the kinetic energy of
these is regarded as zero (equation :eq:`freedom_eq`).

CHARMM Shell Model Self-Induction
---------------------------------

.. index:: single: WWW
   
The `CHARMM <https://www.charmm.org/>`__ model for self-induced
polarisablility relies on the Druder formalism as described above.
However, the CHARMM core-shell model interactions, conventions and
controls have further specificity :cite:`mackerell-16a` that
is worth outlining in DL_POLY_4 terms.

To enable the CHARMM core-shell model in , the user, *at the very
least*, needs to specify atomic polarisabilities (and, optionally, the
respective Thole dumping factors) for all cores in the MPOLES file (see
Section :ref:`mpoles-file`), make sure that reading
MPOLES is triggered by using the directive in the FIELD file (see
Section :ref:`field-file`), and use the directive in
CONTROL (see Section :ref:`control-file`). 

.. note::

   If no Thole dumping factors are specified in MPOLES as well as no
   global Thole dumping factor is (optional) provided with the above
   directive in CONTROL, then a default one of 1.3 is assumed for all
   inducible particles! Also, if no Thole dumping factors are specified in
   MPOLES and a **zero** global Thole dumping factor is (optional) provided
   in CONTROL that will also invalidate the use of CHARMM scaled
   electrostatics in the simulation although the option will help with
   checking, verifying and setting CHARMM related defaults for shell
   charges, core-shell spring force constants and atomic polarisations!

Equation :eq:`druder_eq` governs the relation between the
force constant, :math:`k_{2}` (positive), the shell charge
:math:`q_{s}`, and the atomic polarisability, :math:`\alpha` (positive).
Thus if one is missing, undefined or zero, it can be recovered from the
rest in . CHARMM only allows for :math:`q_{s}` to be recovered. **Note**
that in DL_POLY_4 if any :math:`q_{s}` is recovered, it has an opposite
sign to that of its corresponding :math:`q_{c}`! In the special case
when all Druders’ force constants, :math:`k_{2}`, are undefined or zero
(in FIELD), DL_POLY_4 will resort to using the value of
1000 kcal mol\ :math:`^{-1}`\ Å\ :math:`^{-2}` for all core-shell units,
as per CHARMM recommendations in :cite:`mackerell-16a`. In
all other cases if two (or more) of these three quantities are undefined
or zero DL_POLY_4 will terminate execution in a controlled manner,
indicating the exact nature of the problem.

CHARMM postulates a scaled intra-molecular self-induction between 1-2
and 1-3 intra-molecular neighbours. Thus cross core-shell coulombic
interactions within *conventional* bonded interactions are not fully
excluded (dipole-dipole only, no charge-dipole). They are scaled (see
below) for all possible 1-2 (chemical bonds or constraints) and 1-2-3
(chemical bond angles) qualifying neighbours. For example, let 1-2-3-4
define a torsion angle, then all 1-2-3 and 2-3-4 core-shell
cross-interaction are considered. 

.. note:: 
   **Note** that this is not the case for
   the 1-4 one, which is excluded in this scenario but it could still be
   scaled via torsion 1-4 coulombic scaling!. In DL_POLY_4’s CHARMM context, also,
   then all possible coulombic cross-interactions are considered in the
   case of an inversion angle interaction. It is worth **stressing** that
   in DL_POLY_4 will further exclude (disregard as non-contributing) any
   core-core bonded interactions if they are frozen or mapped onto a RB!

Let’s demonstrate this for a conventional bond angle unit, 1-2-3, with
only members 1 and 2 being polarisable (having a core and a shell). So,
the only CHARMM scaled intra-molecular coulombic interactions for this
scenario will be the four pairs :math:`1_{core}`-\ :math:`2_{core}`,
:math:`1_{core}`-\ :math:`2_{shell}`, :math:`1_{shell}`-\ :math:`2_{core}`
and :math:`1_{shell}`-\ :math:`2_{shell}` (dipole-dipole interactions
only, the charge-dipole interactions are not considered!). In the case
when members 1 and 3 are frozen or in a RB configuration then the pair
:math:`1_{core}`-\ :math:`2_{core}` will be disregarded from that list.

The CHARMM, self-induced intra-molecular (dipole-dipole) interactions
are scaled by a factor of

.. math::

   S(r_{ij}) = 1 -\left[ 1 + \frac{1}{2}
   \frac{a_{i}+a_{j}}{(\alpha_{i}\alpha_{j})^{1/6}} r_{ij}\right]~
   \exp\left(-\frac{a_{i}+a_{j}}{(\alpha_{i}\alpha_{j})^{1/6}} r_{ij}\right)~~,

where :math:`r_{ij}` is the distance between atoms i and j,
:math:`\alpha_{i}` and :math:`\alpha_{j}` are the respective atomic
polarisabilities, and :math:`a_{i}` and :math:`a_{j}` are the respective
atomic (Thole) damping constants. This equation is equivalent to a
smeared charge distribution described by

.. math:: \rho = \frac{a^{3}}{8\pi}~\exp\left(-\frac{a}{(\alpha_{i}\alpha_{j})^{1/6}} r_{ij}\right)~~,

with :math:`a~=~a_{i}+a_{j}`, as originally proposed by Thole
:cite:`thole-81a`.

Dynamical (Adiabatic Shells) Shell Model
----------------------------------------

The dynamical shell model is a method of incorporating polarisability
into a molecular dynamics simulation. The method used in DL_POLY_4 is
that devised by Fincham *et al* :cite:`fincham-93a` and is
known as the adiabatic shell model.

In the adiabatic method, a fraction of the atomic mass is assigned to
the shell to permit a dynamical description. The fraction of mass,
:math:`x`, is chosen to ensure that the natural frequency of vibration
:math:`\nu_{\tt core-shell}` of the harmonic spring (which depends on
the reduced mass, i.e.

.. math:: \nu_{\tt core-shell} = \frac{1}{2\pi} \left[ \frac{k_{2}}{x(1-x)m} \right]^{1/2}~,

with :math:`m` the rigid ion atomic mass) is well above the frequency of
vibration of the whole atom in the bulk system. Dynamically, the
core-shell unit resembles a diatomic molecule with a harmonic 
:index:`bond<potential;chemical bond>`,
however, the high vibrational frequency of the bond prevents effective
exchange of kinetic energy between the core-shell unit and the remaining
system. Therefore, from an initial condition in which the core-shell
units have negligible internal vibrational energy, the units will remain
close to this condition throughout the simulation. This is essential if
the core-shell unit is to maintain a net :index:`polarisation`. (In practice,
there is a slow leakage of kinetic energy into the core-shell units, but
this should should not amount to more than a few percent of the total
kinetic energy. To determine safe shell masses in practice, first a
rigid ion simulation is performed in order to gather the velocity
autocorrelation functions, VAF, of the ions of interest to polarise.
Each VAF is then Fast Fourier transformed to find their highest
frequency of interaction, :math:`\nu_{\tt rigid-ion}`. It is, then, a
safe choice to assign a shell mass, :math:`x~m`, so that
:math:`\nu_{\tt core-shell} \geq 3~\nu_{\tt rigid-ion}`. The user must
make sure to assign the correct mass, :math:`(1-x)~m`, to the core!)

Relaxed (Massless Shells) Model
-------------------------------

The relaxed shell model is presented in :cite:`lindan-93a`,
where shells have no mass and as such their motion is not governed by
the usual Newtonian equation, whereas their cores’ motion is. Because of
that, shells respond instantaneously to the motion of the cores: for any
set of core positions, the positions of the shells are such that the
force on every shell is zero. The energy is thus a minimum with respect
to the shell positions. This represents the physical fact that the
system is always in the ground state with respect to the electronic
degrees of freedom.

Relaxation of the shells is carried out at each time step and involves a
search in the multidimensional space of shell configurations. The search
in DL_POLY_4 is based on the powerful conjugate-gradients technique
:cite:`shewchuk-94a` in an adaptation as shown in
:cite:`lindan-93a`. Each time step a few iterations
(10:math:`\div`\ 30) are needed to achieve convergence to zero net
force.

Breathing Shell Model Extension
-------------------------------

While for low-symmetry structures, the conventional, dipolar rigid shell
model (RSM) is sufficient to absorb most of the effects of partial
covalency/ionic polarisation, for some high-symmetry systems, a
breathing shell model (BSM) :cite:`schroder-66a` is used as
a refinement to represent the contribution of higher-order charge
deformations (of oxide species). This is done by the inclusion of
non-central ion interaction to account for a finite ion shell radius,
:math:`r_{i}^{o}`, which is allowed to deform isotropically under its
environment. However, all short-range repulsion potentials, i.e.
**vdw**, a BSM ion interacts by with its environment, must act upon the
radius of the ion, :math:`U=U(r_{i}-r_{i}^{o})`, rather than the nuclear
position, :math:`U=U(r_{i})`! A further constraining potential is then
added to represent the self-energy of the ion’s breathing shell. This
most commonly uses the same shape as the one of the harmonic bond, see
Section :ref:`bond-potentials`:

.. math:: U(r_{ij}) = \frac{1}{2} k (r_{ij}^{BSM}-r_{ij}^{o})^{2}~~,

where :math:`i` and :math:`j` are the intramolecular index of the ion’s
core and shell respectively. Hence, to employ the BSM, the user needs to
specify an extra bonded interaction for each BSM’ed core-shell pair in
the relevant **bonds** sections in the FIELD file for all molecules that
contain BSM ions. It is worth noting that the BSM energy and virial are
thus part of the bonds’ energy and virial and their calculation as part
of the bonds forces routine bonds_forces.

The most significant consequence of the introduction of the BSM is that,
by coupling the repulsive interactions via common shell radii, it
creates a many-body effect that is able to reproduce the Cauchy
violation (:math:`C_{44} \ne C_{12}`) for rock salt structured
materials.

Further Notes
-------------

In DL_POLY_4 the core-shell forces of the rigid shell model are handled
by the routine core_shell_forces. In case of the adiabatic shell model
the kinetic energy is calculated by core_shell_kinetic and temperature
scaling applied by routine core_shell_quench. In case of the relaxed
shell model shell are relaxed to zero force by core_shell_relaxed.

.. note::
   
   DL_POLY_4 determines which shell model to use by scanning
   shell weights provided the FIELD file (see
   Section :ref:`field-file`). If all shells have zero
   weight the DL_POLY_4 will choose the relaxed shell model. If no shell
   has zero weight then DL_POLY_4 will choose the dynamical one. In case
   when some shells are massless and some are not DL_POLY_4 will terminate
   execution controllably and provide information about the error and
   possible possible choices of action in the OUTPUT file (see
   Section :ref:`output-file`).

.. note::
   
   All DL_POLY_4’s :index:`shell models<polarisation;shell models>` 
   can be used in conjunction with the
   methods for long-ranged forces described above. This also includes uses
   in the context of multipolar electrostatics where self-induced
   polarisation is often a crucial part of polarisable force-field models
   such as CHARMM, AMBER, AMOEBA. 
   
Currently, there are the following
restrictions within DL_POLY_4:

-  shell particles are restricted to only bear a charge.

-  shell particles cannot be frozen, part of a rigid body formation,
   constraint bonded, PMF constrained or tethered.

-  shell particles cannot have shells.

-  a core and shell unit cannot be in a relation via an angle, dihedral
   or inversion type of interaction. However, they can be in a bond type
   of interactions (see the BSM section above).

It is worth noting that bonded interactions such as chemical bonds,
angles, dihedrals and inversions can be defined over any possible
mixture of core and shell members of different core-shell units! For
example, in the solid state materials community it is common to define
bonds and angles over the shells as the shells represent the electronic
clouds between which the bonding occurs. However, this is right the
opposite in the liquid, organic and bio-chemical communities, where all
bonding is between the nuclei (cores) and the shells are only there to
account purely for the polarisability. As DL_POLY_4 allows for too much
flexibility in the space of possible model definitions, it is strongly
advised that the modeller must exercise great care to define the correct
force-field model representation within DL_POLY_4 input files!