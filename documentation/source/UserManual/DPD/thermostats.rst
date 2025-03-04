.. _DPD_thermostats:

Thermostats
===========

DL_POLY_5 supports several DPD based thermostats which can be applied to 
both DPD or traditional MD type simulations. To use any DPD based thermostat the ``ensemble_method``
directive in the ``CONTROL`` file must be set to ``dpd``. The ``ensemble_dpd_order`` directive is then 
used to specify the type of DPD based thermostat to be used. 

All DPD based thermostats require setting the value of the drag coefficient :math:`\gamma_{ij}` 
from :eq:`DPD_drag_force_eq`. This can be specified in the ``FIELD`` file as a per-pair interaction 
parameter when using a DPD based force field.
If the drag coefficient has not been specified in the ``FIELD`` file as part of a dpd based interaction, 
the ``ensemble_dpd_drag`` directive can be used to set a global value for the drag coefficient to be used 
within the DPD thermostats.

.. note::
    If the ``ensemble_dpd_drag`` value is specified then all interactions will use this value regardless 
    of any values specified in the ``FIELD`` file. 

Molecular Dynamics Velocity Verlet for DPD
------------------------------------------

Here a standard Velocity Verlet integration scheme is employed with the random and drag forces, :math:`\underline{F}_{ij}^{\text{R}}` 
and :math:`\underline{F}_{ij}^{\text{D}}` from :eq:`DPD_eq` used as the system thermostat, i.e. the thermostatting force 
:math:`\underline{F}_{ij}^{\text{T}} = \underline{F}_{ij}^{\text{D}} + \underline{F}_{ij}^{\text{R}}`. This thermostatting force is 
combined with all other forces between particles; pairwise conservative (standard and/or density-dependent), bonding, electrostatic, 
external (body) forces etc. 

The combination of the DPD thermostat with the standard molecular dynamics type Velocity Verlet algorithm is the simplest and least 
time-consuming DPD based thermostatting algorithm available in DL_POLY_5. The drag force does, however, depend upon particle velocities
and is therefore only approximated using the mid-step values. This can produce a system temperature that is higher than that specified by 
the user and would require a small time step :math:`\Delta t` to reduce the offset to a tolerable level.  

To use this thermostat, in the ``CONTROL`` file, set the ``ensemble_method`` and ``ensemble_dpd_order`` directives to ``dpd`` and ``mdvv`` 
respectively.  

Shardlow Splitting
------------------

Shardlow splitting :cite:`shardlow-03a` exploits an operator splitting approach to integrate 
the drag and random, forces from :eq:`DPDF_eq`, in the DPD thermostat. The changes in 
particle velocities due to the thermostat during a timestep can be determined separately
from changes due to other forces (conservative, bonding interactions etc.) by rigorously 
expanding velocity Verlet integration of the thermostatting forces. The separate integration 
of drag and random forces can be described as a Trotter (first-order) or Strang (second-order) 
operator. Both first-order and second-order Shardlow operators can be coupled with a barostat to 
allow for constant pressure, surface area or surface tension simulations :cite:`lisal2011`. 

First-order 
~~~~~~~~~~~

The **first-order** Shardlow operator can be achieved in two stages, equivalent to the first and second 
stages of velocity Verlet integration of drag and random forces (:math:`\underline{f}_{ij}^{D}` and 
:math:`\underline{f}_{ij}^{R}`). The first stage is simply integration of the forces over half of the 
timestep :math:`\Delta t` for all particle pairs(:math:`i` and :math:`j`) within the thermostat cutoff:

.. math:: 
    :label: shardlow_stage1_eq 

    \underline{v}_i (t + \tfrac{1}{2}\Delta t) &= \underline{v}_i (t) - \frac{\Delta t}{2 m_i} \left[ - \gamma_{ij} w^D \left( r_{ij} \right)\left(\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij} (t) \right) + \sigma_{ij} w^R \left(r_{ij} \right) \zeta_{ij} \Delta t^{-\frac{1}{2}} \right] \frac{\underline{r}_{ij}}{r_{ij}} \\
    \underline{v}_j (t + \tfrac{1}{2}\Delta t) &= \underline{v}_j (t) + \frac{\Delta t}{2 m_j} \left[ - \gamma_{ij} w^D \left( r_{ij} \right)\left(\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij} (t) \right) + \sigma_{ij} w^R \left(r_{ij} \right) \zeta_{ij} \Delta t^{-\frac{1}{2}} \right] \frac{\underline{r}_{ij}}{r_{ij}}.

The equations for the second velocity Verlet stage would ordinarily be implicit due to the inclusion of the
relative velocity between particle pairs at the end of the timestep in the dissipative force. However, these 
equations can be rewritten in an explicit form that gives the velocity changes based on the mid-step relative 
velocities:

.. math:: 
    :label: shardlow_stage2_eq

    \underline{v}_i (t + \Delta t) &= \underline{v}_i (t + \tfrac{1}{2}\Delta t) - \frac{\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}}{2 m_{i}}\frac{\underline{r}_{ij}}{r_{ij}} \\
    &+ \frac{\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t}{2m_{i}\left(1 + \frac{1}{2}\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t\mu_{ij}^{-1}\right)} \left[\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij}\left(t + \frac{1}{2}\Delta t\right) + \frac{1}{4}\mu_{ij}^{-1}\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}\right]\frac{\underline{r}_{ij}}{r_{ij}} \\
    \underline{v}_j (t + \Delta t) &= \underline{v}_j (t + \tfrac{1}{2}\Delta t) - \frac{\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}}{2 m_{i}}\frac{\underline{r}_{ij}}{r_{ij}} \\
    &+ \frac{\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t}{2m_{i}\left(1 + \frac{1}{2}\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t\mu_{ij}^{-1}\right)} \left[\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij}\left(t + \frac{1}{2}\Delta t\right) + \frac{1}{4}\mu_{ij}^{-1}\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}\right]\frac{\underline{r}_{ij}}{r_{ij}}

where :math:`\mu_{ij} = \frac{m_i m_j}{m_i + m_j}` is the reduced mass between the two particles, i.e. 
:math:`\frac{1}{2}\mu_{ij}^{-1}` is their mean reciprocal mass. 
Other forces can subsequently be integrated using the standard velocity Verlet algorithm, which do not depend upon the particle velocities. 

First-order shardlow splitting is specified by setting the ``ensemble_dpd_order`` directive to ``first`` or 
``1`` in the ``CONTROL`` file. 

Second-order 
~~~~~~~~~~~~

The **second-order** Shardlow operator is virtually identical to the first-order operator, except for replacing all instances of :math:`\Delta t` with :math:`\frac{\Delta t}{2}` in the above velocity corrections and applying them twice, both before (starting at :math:`t`) and after (starting at :math:`t + \tfrac{1}{2}\Delta t`) integration of other forces. In both cases, correct integration of drag and random forces from known velocity data allows larger timesteps to be used than for standard Velocity Verlet integration. 

Second-order Shardlow splitting is specified by setting the ``ensemble_dpd_order`` directive to ``second`` or 
``1`` in the ``CONTROL`` file. 

Zeroth-order
~~~~~~~~~~~~

An alternative approach is **zeroth-order** Shardlow splitting :cite:`lisal2011` where the typical Verlocity 
Verlet style integration is replaced by an implicit algorithm over a single step, 

.. math::
    :label: zeroth_shardlow_eq 

    \underline{v}_i (t + \Delta t) &= \underline{v}_i (t) - \frac{\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}}{m_{i}}\frac{\underline{r}_{ij}}{r_{ij}} \\
    &+ \frac{\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t}{m_{i}\left(1 + \gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t\mu_{ij}^{-1}\right)} \left[\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij}\left(t \right) + \frac{1}{2}\mu_{ij}^{-1}\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}\right]\frac{\underline{r}_{ij}}{r_{ij}} \\
    \underline{v}_j (t + \Delta t) &= \underline{v}_j (t) - \frac{\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}}{m_{i}}\frac{\underline{r}_{ij}}{r_{ij}} \\
    &+ \frac{\gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t}{m_{i}\left(1 + \gamma_{ij}w^{D}\left(r_{ij}\right)\Delta t\mu_{ij}^{-1}\right)} \left[\frac{\underline{r}_{ij}}{r_{ij}} \cdot \underline{v}_{ij}\left(t \right) + \frac{1}{2}\mu_{ij}^{-1}\sigma_{ij}w^{R}\left(r_{ij}\right)\zeta_{ij}\Delta t^{\frac{1}{2}}\right]\frac{\underline{r}_{ij}}{r_{ij}}


Second-order Shardlow splitting is specified by setting the ``ensemble_dpd_order`` directive to ``zeroth`` or 
``0`` in the ``CONTROL`` file. 
