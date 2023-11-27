Introduction
============

As a default the DL_POLY_4integration :index:`algorithm`s are based on the
Velocity :index:`Verlet<algorithm;Verlet>` (VV) scheme, which is both simple and time reversible
:cite:`allen-89a`. It generates trajectories in the
microcanonical (NVE) ensemble in which the total energy (kinetic plus
potential) is conserved. If this property drifts or fluctuates
excessively in the course of a simulation it indicates that the timestep
is too large or the potential cutoffs too small (relative r.m.s.
fluctuations in the total energy of :math:`10^{-5}` are typical with
this algorithm).

The VV algorithm has two stages (VV1 and VV2). At the first stage it
requires values of position (:math:`\underline{r}`), velocity
(:math:`\underline{v}`) and force (:math:`\underline{f}`) at time :math:`t`. The
first stage is to advance the velocities to :math:`t+(1/2)\Delta t` by
integration of the force and then to advance the positions to a full
step :math:`t+\Delta t` using the new half-step velocities:

#. VV1:

   .. math::

      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow \underline{v}(t) + {\Delta t \over 2} \;
      {\underline{f}(t) \over m}~~,

   where :math:`m` is the mass of a site and :math:`\Delta t` is the
   timestep

   .. math::

      \underline{r}(t + \Delta t) \leftarrow \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2}\Delta t)

#. FF: Between the first and the second stage a recalculation of the
   force at time :math:`t+\Delta t` is required since the positions have
   changed

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2: In the second stage the half-step velocities are advanced to to
   a full step using the new force

   .. math::

      \underline{v}(t + \Delta t) \leftarrow  \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; {\underline{f}(t + \Delta t) \over m}

The instantaneous kinetic energy, for example, can then be obtained from
the atomic velocities as

.. math:: E_{kin}(t) = {1 \over 2} \sum_{1}^{\cal N} m_{i} v_{i}^{2}(t)~~, \label{ekin}
   :label: ekin_eq

and assuming the system has no net momentum the instantaneous
temperature is

.. math:: {\cal T}(t) = \frac{2}{k_{B} f} E_{kin}(t)~~, \label{tinst}
   :label: tinst_eq

where :math:`i` labels particles (that can be free atoms or :index:`rigid<rigid body>`
bodies), :math:`{\cal N}` the number of particles (free atoms and rigid
bodies) in the system, :math:`k_{B}` the Boltzmann’s constant and
:math:`f` the number of degrees of freedom in the system.

.. math::
   :label: freedom_eq

   f = 3{\cal N} - 3{\cal N}_{frozen} - 3{\cal N}_{shells} -
   {\cal N}_{constraints} - 3 - p~~. \label{freedom}

Here :math:`{\cal N}_{frozen}` indicates the number of frozen atoms in
the system, :math:`{\cal N}_{shells}` number of core-shell units and
:math:`{\cal N}_{constraints}` number of bond and PMF constraints. Three
degrees of freedom are subtracted for the centre of mass zero net
momentum (which we impose) and :math:`p` is zero for periodic or three
for non-periodic systems, where it accounts for fixing angular momentum
about origin (which we impose).

In the case of :index:`rigid<rigid body>` bodies (see Section \ :ref:`rigid`) the first
part of equation :eq:`freedom_eq`

.. math:: f\prime = 3{\cal N} - 3{\cal N}_{frozen}

splits into

.. math::

   f\prime = \left( 3{\cal N}^{FP} - 3{\cal N}^{FP}_{frozen} \right) +
             \left( 3{\cal N}^{RB \mathbf{ (tra)}} - 3{\cal N}^{RB \mathbf{ (tra)}}_{frozen} \right) +
             \left( 3{\cal N}^{RB \mathbf{ (rot)}} - 3{\cal N}^{RB \mathbf{ (rot)}}_{frozen} \right)

or

.. math:: f\prime = f^{FP} + f^{RB \mathbf{ (tra)}} + f^{RB \mathbf{ (rot)}}~~.

Here FP stands for a free particle, i.e. a particle not participating in
the constitution of a rigid body, and RB for a rigid body. In general a
rigid body has 3 translational (:math:`\mathbf{ tra}`) degrees of
freedom, corresponding to its centre of mass being allowed to move in
the 3 general direction of space, and 3 rotational
(:math:`\mathbf{ rot}`), corresponding to the RB being allowed to rotate
around the 3 general axis in space. It is not far removed to see that
for a not fully frozen rigid body one must assign 0 translational
degrees of freedom but depending on the "frozenness" of the RB one may
assign 1 rotational degrees of freedom when all the frozen sites are in
line (i.e. rotation around one axis only) or 3 when just one site is
frozen.

The routine ``nve_0_vv`` implement the Verlet algorithm in velocity
:index:`verlet<algorithm;Verlet>` for free particles and calculate the instantaneous temperature.
Whereas the routines ``nve_1_vv`` implements the same for systems also
containing rigid bodies. The conserved quantity is the total energy of
the system

.. math:: {\cal H}_{\rm NVE} = U + E_{kin}~~,

where :math:`U` is the potential energy of the system and
:math:`E_{kin}` the kinetic energy at time :math:`t`.

The full selection of integration algorithms within DL_POLY_4is as
follows:

.. list-table::

   * - ``nve_0_vv`` 
     - Constant E :index:`algorithm<ensemble;NVE>`
   * - ``nve_1_vv``
     - The same as the above but also incorporating RB integration
   * - ``dpd_thermostat``
     - Constant T :index:`algorithm<ensemble;DPD NVT>` (DPD:cite:`shardlow-03a`)
   * - ``nvt_e0_vv`` 
     - Constant :math:`E_{kin}` :index:`algorithm<ensemble;Evans NVT>` (Evans:cite:`evans-84a`)
   * - ``nvt_e1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``nvt_l0_vv`` 
     - Constant T :index:`algorithm<ensemble;Langevin NVT>` (Langevin:cite:`adelman-76a`)
   * - ``nvt_l1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``nvt_l2_vv`` 
     - Constant T :index:`algorithm<ensemble;Inhomogeneous Langevin NVT>` (inhomogeneous Langevin :cite:`duffy-07a`)
   * - ``nvt_a0_vv`` 
     - Constant T :index:`algorithm<ensemble;Anderson NVT>` (Andersen :cite:`andersen-79a`)
   * - ``nvt_a1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``nvt_b0_vv`` 
     - Constant T :index:`algorithm<ensemble;Berendsen NVT>` (Berendsen :cite:`berendsen-84a`)
   * - ``nvt_b1_vv``
     - The same as the above but also incorporating RB integration
   * - ``nvt_h0_vv`` 
     - Constant T :index:`algorithm<ensemble;Nosé-Hoover NVT>` (Hoover :cite:`hoover-85a`)
   * - ``nvt_h1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``nvt_g0_vv`` 
     - Constant T :index:`algorithm<ensemble;Gentle Stochastic NVT>` (GST :cite:`leimkuhler-09a`)
   * - ``nvt_g1_vv``
     - The same as the above but also incorporating RB integration
   * - ``npt_l0_vv`` 
     - Constant T,P :index:`algorithm<ensemble;Langevin NPT>` (Langevin :cite:`quigley-04a`)
   * - ``npt_l1_vv``
     - The same as the above but also incorporating RB integration 
   * - ``npt_b0_vv`` 
     - Constant T,P :index:`algorithm<ensemble;Berendsen NPT>` (Berendsen :cite:`berendsen-84a`)
   * - ``npt_b1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``npt_h0_vv`` 
     - Constant T,P :index:`algorithm<ensemble;Nosé-Hoover NVT>` (Hoover :cite:`hoover-85a`)
   * - ``npt_h1_vv``
     - The same as the above but also incorporating RB integration
   * - ``npt_m0_vv`` 
     - Constant T,P :index:`algorithm<ensemble;Martyna-Tuckerman-Klein NPT>` (Martyna-Tuckerman-Klein :cite:`martyna-96a`)
   * - ``npt_m1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``npt_l0_vv``
     -  Constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Langevin N$\simga$T>` (Langevin:cite:`quigley-04a`)
   * - ``npt_l1_vv``
     - The same as the above but also incorporating RB integration
   * - ``nst_b0_vv`` 
     - Constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Berendsen N$\sigma$T>` (Berendsen :cite:`berendsen-84a`)
   * - ``nst_b1_vv`` 
     - The same as the above but also incorporating RB integration
   * - ``nst_h0_vv`` 
     - Constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Nosé-Hoover N$\sigma$T>` (Hoover :cite:`hoover-85a`)
   * - ``nst_h1_vv`` 
     - The same as the above but also incorporating RBs integration
   * - ``nst_m0_vv`` 
     - Constant T,\ :math:`\underline{\underline{\mathbf{\sigma}}}` :index:`algorithm<ensemble;Martyna-Tuckerman-Klein N$\sigma$T>` (Martyna-Tuckerman-Klein :cite:`martyna-96a`)
   * - ``nst_m0_vv`` 
     - The same as the above but also incorporating RB integration

It is worth noting that the last four ensembles are also optionally
available in an extended from to constant normal pressure and constant
surface area, NP\ :math:`_{n}`\ AT, or constant surface tension,
NP\ :math:`_{n}\gamma`\ T :cite:`ikeguchi-04a`.
