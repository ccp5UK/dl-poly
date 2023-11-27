.. _DPD-all:

Appendix A: DL_POLY_4 Dissipative Particle Dynamics
+++++++++++++++++++++++++++++++++++++++++++++++++++

Introduction
============

Although in a Molecular Dynamics sense Dissipative Particle Dynamics
(DPD) is regarded as a type of thermostat, on its own it is an
off-lattice, discrete particle method for modelling mesoscopic systems
in a fluid state. It has little in common with Lattice Gas Automata and
Lattice Boltzmann dynamics methods :cite:`hoogerbrugge-92a`,
except in its application to systems of similar length and time scales,
and last but not least that it captures hydrodynamics behaviour.

The DPD method inherits its methodology from Brownian/Langevin Dynamics
(BD). However, it differs from BD in an important way: it is *Galilean
invariant* and for this reason conserves hydrodynamic behaviour, while
the BD method does not (its microscopic behaviour is only diffusive.
Many systems in their fluid state are crucially dependent on
hydrodynamic interactions and it is essential to retain this feature in
their models. DPD is particularly useful for simulating coarse-grained
systems on the near-molecular scale, such as polymers, biopolymers,
lipids, emulsions and surfactants – systems in which large scale
structure evolves on a time scale that is too long to be modelled
effectively by traditional MD. The particles in DPD
:cite:`espanol-95a` are not regarded as molecules in a fluid
but as lumps of molecules grouped to form a *fluid particle* in much the
same spirit as the renormalisation group has been applied in polymer
physics where lumps of monomers are grouped to form a *bead*. Hence, the
beads are regarded as carriers of momentum.

It is worth noting that DPD may also be used when such systems
experience shear and flow gradients.

Outline of Method
=================

Following :cite:`groot-97a` the DPD algorithm can be
summarised by the following:

-  A condensed phase system may be modelled as a system of ‘free’
   particles interacting directly through *soft* forces. Note that
   DL_POLY_4 allows for the application of the DPD thermostat beyond
   systems of free particles only. Thus it will be valid on systems with
   any inratamolecular like interactions.

-  The system is coupled to a heat bath via stochastic forces, which act
   on the particles in a pairwise manner.

-  The particles also experience a damping or drag force, which also
   acts in a pairwise manner.

-  Thermodynamic equilibrium is maintained through the balance of the
   stochastic and drag forces, i.e. the method satisfies the
   fluctuation-dissipation theorem.

-  At equilibrium (or steady state) the properties of the system are
   calculated as averages over the individual particles, as in
   traditional Molecular Dynamics.

Therefore, the equation of motion are the same as these for the
microcanonical ensemble (NVE) but force, :math:`f_{i}`, on particle
:math:`i` is now a sum of pair forces:

.. math:: \underline{f}_i = \sum_{j \neq i}^N \left( \underline{f}_{ij}^{C} + \underline{f}_{ij}^{D} + \underline{f}_{ij}^{R} \right)~~, \label{DPD}
   :label: DPD_eq

in which :math:`\underline{f}_{ij}^{C}`, :math:`\underline{f}_{ij}^{D}`
and :math:`\underline{f}_{ij}^{R}` are the *conservative*, *drag* and
*random* (or *stochastic*) pair forces respectively. Each represents the
force exerted on particle :math:`i` due to the presence of particle
:math:`j`.

The conservative interactions are usually *soft* (i.e. weakly
interacting) so that the particles can pass by each other (or even
through each other) relatively easily so that equilibrium is achieved
quickly. A common form of interaction potential is an inverse parabola
(part of VDW types of potentials, see Section \ :ref:`vdw`)

.. math::
   :label: DPDU_eq

   V(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{A_{ij}}{2}~r_{c}~\left(1-\frac{r_{ij}}{r_{c}}\right)^{2} & r_{ij} < r_{c} \\
   0 & r_{ij} \ge r_{c} \end{array} \right.~~, \label{DPDU}

where :math:`r_{ij} = |\underline{r}_{j}-\underline{r}_{i}|`,
:math:`r_{c}` is a cutoff radius and :math:`A_{ij}` is the interaction
strength (that may be the same for all particle pairs or may be
different for different particle types).

Equation :eq:`DPD_eq` gives rise to a repulsive force of the
form:

.. math::
   :label: DPDF_eq

   \underline{f}_{ij}^{C} = A_{ij}~w^{C}(r_{ij}) \frac{\underline{r}_{ij}}{r_{ij}} =
   A_{ij} \left( 1 - \frac{r_{ij}}{r_{c}} \right) \frac{\underline{r}_{ij}}{r_{ij}}~~. \label{DPDF}

This is the deterministic or *conservative* force
:math:`\underline{f}_{ij}^{C}` exerted on particle :math:`i` by particle
:math:`j`. Note the switching function:

.. math::
   :label: DPDS_eq

   w^{C}(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
   \left(1-\frac{r_{ij}}{r_{c}}\right) & r_{ij} < r_{c} \\
   0 & r_{ij} \ge r_{c} \end{array} \right.~~, \label{DPDS}

and the force are zero when :math:`r_{ij} \ge r_{c}` and thus the
particles have an effective diameter of :math:`1` in units of the cutoff
radius :math:`r_{c}`. In the DL_POLY_4 context all inter- and
intra-molecular forces will fall into this category of force!

The stochastic forces experienced by the particles is again pairwise in
nature and takes the form:

.. math:: \underline{f}_{ij}^{R} = \sigma_{ij} w^{R}(r_{ij}) \zeta_{ij} \Delta t^{-\frac{1}{2}} \frac{\underline{r}_{ij}}{r_{ij}}~~,

in which :math:`\Delta t` is the time step and :math:`w^{R}(r_{ij})` is
a switching function which imposes a finite limit on the range of the
stochastic force. :math:`\zeta_{ij}` is a random number with zero mean
and unit variance. The constant :math:`\sigma_{ij}` is related to the
temperature, as is understood from the role of the stochastic force in
representing a heat bath.

Finally, the particles are subject to a drag force, which depends on the
relative velocity between interacting pairs of particles:

.. math::

   \underline{f}_{ij}^{D} = -\gamma_{ij} w^{D}(r_{ij})
   \left(\underline{r}_{ij} \cdot \underline{v}_{ij}\right) \frac{\underline{r}_{ij}}{r_{ij}^2}~~,

where :math:`w^{D}(r_{ij})` is once again a switching function and
:math:`\underline{v}_{ij} = \underline{v}_{j}-\underline{v}_{i}` is the
inter-particle relative velocity. The constant :math:`\gamma_{ij}` is
the drag coefficient. It follows from the fluctuation-dissipation
theorem that for thermodynamic equilibrium to result from this method
the following relations must hold:

.. math:: \sigma_{ij}^2 = 2~\gamma_{ij} k_{B} T \label{DPDC1} 
   :label: DPDC1 

.. math:: w^{D}(r_{ij}) = \left[w^{R}(r_{ij})\right]^{2}. \label{DPDC2} 
   :label: DPDC2

In practice, the switching functions are defined through:

.. math:: w^{R}(r_{ij}) = \left[w^{C}(r_{ij})\right]^{2}~~,

which ensures that all interactions are switched off at the range
:math:`r_{ij} = r_{c}`.

In many DPD simulations, the stochastic and drag coefficients are often
constant for all interactions, i.e. :math:`\sigma_{ij} \equiv \sigma`
and :math:`\gamma_{ij} \equiv \gamma`, although this assumption does not
have to apply. In DL_POLY_4 the :math:`\gamma_{ij}` coefficients may be
supplied at the end of each specified vdw interaction potential as a
parameter further to the last one for the particular vdw potential form.
For a DPD thermostat to work correctly all possible two body
interactions must be defined and all :math:`\gamma_{ij} \ne 0`. What
DL_POLY_4 will attempt first, if a two body interaction is missing, is
to derive it using mixing rules (default may be overridden by user
specification). However, if any of :math:`\gamma_{ij} = 0` then
DL_POLY_4 will check for the existence of a global :math:`\gamma` that
may be optionally supplied by the user on the **ensemble nvt dpd**\ INT
line and if it is non-zero a global override will occur. Otherwise, when
the requirements for a DPD thermostat are not satisfied, everything else
will result in a controlled termination.

Equation of state and dynamic properties
========================================

The form of the conservative force determines the equation of state for
a DPD fluid, which can be derived using the virial theorem to express
system pressure as follows:

.. math::

   \begin{aligned}
   {\cal P} =& \rho k_{B}T + \frac{1}{3V} \left\langle \sum_{j>i} (\underline{r}_i-\underline{r}_j) \cdot \underline{f}_{ij}^{C} \right\rangle \\
   =& \rho k_{B} T + \frac{2 \pi}{3} \rho^{2} \int_{0}^{r_{c}} A \left( 1 - \frac{r}{r_{c}} \right) r^{3} g(r)~dr~~,
   \end{aligned}

where :math:`g(r)` is a radial distribution function for the soft sphere
model :cite:`groot-97a` and :math:`\rho` is the DPD particle
density. For sufficiently large densities (:math:`\rho > 2`),
:math:`g(r)` takes the same form and the equation of state can be
well-approximated by:

.. math:: {\cal P} = \rho k_{B}T + \alpha A \rho^{2}~~,

where the parameter :math:`\alpha \approx 0.101 \pm 0.001` has units
equivalent to :math:`r_{c}^{4}`. This expression permits the use of
fluid compressibilities to obtain conservative force parameters for bulk
fluids, e.g. for water :math:`A \approx 75 k_{B} T/\rho`. Alternative
equations of state may be obtained by modifying the functional form of
conservative interactions to include localized densities (i.e. many-body
DPD) :cite:`pagonabarraga-01a,trofimov-02a`.

Transport coefficients for a DPD fluid without conservative forces can
be derived using the expressions for the drag and stochastic
forces:cite:`groot-97a,koelman-93a,marsh-97a`. The kinematic
viscosity can be found to be

.. math:: \nu \approx \frac{45 k_{B} T}{4 \pi \gamma \rho r_{c}^{3}} + \frac{2 \pi \gamma \rho r_{c}^{5}}{1575}~~,

while the self-diffusion coefficient is given as

.. math:: D \approx \frac{45 k_{B} T}{2 \pi \gamma \rho r_{c}^{3}}.

The ratio of these two properties, the Schmidt number
(:math:`\textnormal{Sc} = \nu / D`), is therefore:

.. math:: \textnormal{Sc} \approx \frac{1}{2} + \frac{(2 \pi \gamma \rho r_{c}^{4})^{2}}{70875 k_{B} T}

and for values of the drag coefficient and density frequently used in
DPD simulations, this value is of the order of unity, which is an
appropriate magnitude for gases but three orders of magnitude too small
for liquids.

This property of standard DPD does *not* rule it out for simulations of
liquid phases except when hydrodynamics are important. It may also be
argued that the self-diffusion of DPD particles might not correspond to
that of individual molecules and thus a Schmidt number of the order
:math:`10^{3}` is unnecessary for modelling liquids
:cite:`peters-04a`. Alternative thermostats are available in
the DL_MESO :cite:`seaton-13a` (`<http://www.ccp5.ac.uk/DL\_MESO/>`_) 
:index:`package<WWW>`, which can model systems
with higher Schmidt numbers :cite:`lowe-99a,stoyanov-05a`.

Derivation of Equilibrium
=========================

The derivation of the DPD algorithm is based on the Fokker-Planck
equation

.. math:: \frac{\partial \rho}{\partial t} = \mathcal{L} \rho \label{FokkerPlanck}
   :label: FokkerPlanck_eq

where :math:`\rho` is the equilibrium distribution function and
:math:`\mathcal{L}` is the evolution operator, which may be split into
*conservative* and *stochastic+dissipative* parts:

.. math:: \mathcal{L} = \mathcal{L}^{C} + \mathcal{L}^{R+D}

with

.. math::
   :label: DPDEvolution_eq

   \begin{aligned}
   \mathcal{L}^{C} =& -\sum_{i=1}^{N} \frac{\underline{p}_{i}}{m_{i}}
   \frac{\partial}{\partial \underline{r}_{i}} - \sum_{i \neq j}^{N}
   \underline{f}_{ij}^{C} \frac{\partial}{\partial \underline{p}_{i}} \\
   \mathcal{L}^{R+D} =& \sum_{i=1}^{N} \hat{e}_{ij} \cdot \frac{\partial}{\partial \underline{p}_{i}}
   \left[ \frac{\sigma^{2}}{2} \left\{w^{R} \left(r_{ij} \right) \right\}^{2} \hat{e}_{ij} \cdot
   \left\{ \frac{\partial}{\partial \underline{p}_{i}} - \frac{\partial}{\partial \underline{p}_{j}} \right\} +
   \gamma w^{D} \left( \hat{e}_{ij} \cdot \underline{v}_{ij} \right) \right]~~, \label {DPDEvolution}\end{aligned}

where :math:`\hat{e}_{ij} = \frac{\underline{r_{ij}}}{r_{ij}}`.

When :math:`\sigma = \gamma = 0` then
equation :eq:`FokkerPlanck_eq` becomes

.. math:: \frac{\partial \rho}{\partial t} = \mathcal{L}^{C} \rho~~,

for which the equilibrium solution is evidently

.. math::

   \rho^{eq} = \frac{1}{Z} \exp \left( \frac{1}{k_{B} T} \left[ \sum_{i=1}^{N}
   \frac{p_{i}^{2}}{2 m_{i}} + \frac{1}{2} \sum_{j \neq i}^{N} \phi (r_{ij}) \right] \right)

which is, of course, the Boltzmann distribution function for an
equilibrium system. Thus it is apparent that for the simulation based on
equation :eq:`FokkerPlanck_eq` to maintain the same
distribution function, the terms in the operator
:math:`\mathcal{L}^{R+D}` of
equation :eq:`DPDEvolution_eq` must sum to zero. It
follows that the conditions given in equations :eq:`DPDC1`
and :eq:`DPDC2` must apply.

Summary of Dissipative Particle Dynamics
========================================

DPD is a simple method that can be viewed as a novel thermostatting
method for molecular dynamics. All that is required is a system of
spherical particles enclosed in a periodic box undergoing time evolution
as a result of the above forces. It should be noted that all computed
interactions are pairwise, which means that the principle of the
conservation of momentum in the system, or *Galilean invariance*, is
preserved. The conservation of momentum is required for the preservation
of hydrodynamic forces. Therefore, the DPD method is an NVT method that
*preserves hydrodynamics*. The presence of hydrodynamics is important in
annealing defects in ordered mesophases :cite:`gonella-97a`.
Thus DPD has an intrinsic advantage over other methods such as
traditional molecular dynamics, dynamic density functional theory (which
are purely *diffusive*!) or Monte Carlo methods, in trying to evolve a
system towards an ordered thermodynamic equilibrium state.
