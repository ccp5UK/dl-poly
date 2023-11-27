
Thermostats
===========

.. index:: 
   single: ensemble;NVE 
   single: ensemble;canonical 
   single: constraints;Gaussian

The system may be coupled to a heat bath to ensure that the average
system temperature is maintained close to the requested temperature,
:math:`T_{\rm ext}`. When this is done the equations of motion are
modified and the system no longer samples the microcanonical ensemble.
Instead trajectories in the canonical (NVT) ensemble, or something close
to it are generated. DL_POLY_4 comes with seven different thermostats:
Evans (Gaussian constraints) :cite:`evans-84a`, Langevin
(both standard :cite:`adelman-76a,izaguirre-01a` and
inhomogeneous :cite:`duffy-07a` variants), Andersen
:cite:`andersen-79a`, Berendsen
:cite:`berendsen-84a`, Nosé-Hoover (N-H)
:cite:`hoover-85a` and the gentle stochastic thermostat
(GST) :cite:`leimkuhler-09a,samoletov-07a` as well as
Dissipative Particle Dynamics (\ :index:`DPD`) method
:cite:`hoogerbrugge-92a,espanol-95a,groot-97a,shardlow-03a`.
Of these, only the Langevin, N-H, GST and DPD algorithms generate *true*
trajectories in the canonical (NVT) ensemble. The rest will produce
properties that typically differ from canonical averages by
:math:`{\cal O}(1/{\cal N})` :cite:`allen-89a` (where
:math:`\cal N` is the number of particles in the system), as the Evans
algorithm generates trajectories in the (NVE\ :math:`_{kin}`) ensemble.

Evans Thermostat (Gaussian Constraints)
---------------------------------------

Kinetic temperature can be made a constant of the equations of motion by
imposing an additional constraint on the system. If one writes the
equations of motion as:

.. math::

   \begin{aligned}
   {d \underline{r}(t) \over d t} =& \underline{v}(t) \nonumber \\
   {d \underline{v}(t) \over d t} =& {\underline{f}(t) \over m} - \chi (t) \;
   \underline{v}(t)~~,\end{aligned}

the kinetic temperature constraint :math:`\chi` can be found as follows:

.. math::
   :label: Evans_eq

   \begin{aligned}
   \frac{d}{dt} {\cal T} \propto\frac{d}{dt} \left(\frac{1}{2}\sum_{i} m_{i} \underline{v}_{i}^{2} \right) =
   \sum_{i} m_{i} \underline{v}_{i} \cdot \frac{d}{dt} \underline{v}_{i} =& 0 \nonumber \\
   \sum_{i} m_{i} \underline{v}_{i}(t) \cdot
   \left\{ \frac{\underline{f}_{i}(t)}{m_{i}} - \chi (t) \; \underline{v}_{i}(t) \right\} =& 0 \label{Evans} \\
   \chi (t) =& \frac {\sum_{i} \underline{v}_{i}(t) \cdot \underline{f}_{i}(t)} {\sum_{i} m_{i} \underline{v}_{i}^{2}(t)}~~, \nonumber\end{aligned}

where :math:`\cal T` is the instantaneous temperature defined in
equation :eq:`tinst_eq`.

The VV implementation of the Evans algorithm is straight forward. The
conventional VV1 and VV2 steps are carried out as before the start of
VV1 and after the end of VV2 there is an application of thermal
constraining. This involves the calculation of :math:`\chi (t)` before
the VV1 stage and :math:`\chi (t+\Delta t)` after the VV2 stage with
consecutive thermalisation on the unthermostated velocities for half a
timestep at each stage in the following manner:

#. Thermostat VV1

   .. math::

      \begin{aligned}
      \chi (t) &\leftarrow& \frac {\sum_{i} \underline{v}_{i}(t) \cdot \underline{f}_{i}(t)} {2~E_{kin}(t)} \nonumber \\
      \underline{v}(t) &\leftarrow& \underline{v}(t)~\exp \left( -\chi (t) {\Delta t \over 2} \right)~~.\end{aligned}

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) &\leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) &\leftarrow& \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t)\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; \left[{\underline{f}(t + \Delta t) \over m} \right]

#. RATTLE_VV2

#. Thermostat VV2

   .. math::

      \begin{aligned}
      \chi (t + \Delta t) &\leftarrow& \frac {\sum_{i} \underline{v}_{i}(t +
      \Delta t) \cdot \underline{f}_{i}(t + \Delta t)} {2~E_{kin}(t + \Delta t)} \nonumber \\
      \underline{v}(t + \Delta t) &\leftarrow& \underline{v}(t + \Delta t)~
      \exp \left( -\chi (t + \Delta t) {\Delta t \over 2} \right)~~.\end{aligned}

The algorithm is self-consistent and requires no iterations.

The conserved quantity by these algorithms is the system kinetic energy.

The VV flavour of the Gaussian constraints algorithm is implemented in
the DL_POLY_4routines ``nvt_e0_vv``. The routine ``nvt_e1_vv`` implement
the same but also incorporate RB dynamics.

Langevin Thermostat
-------------------

The Langevin thermostat works by coupling every particle to a viscous
background and a stochastic heath bath (Brownian dynamics) such that

.. math::

   \begin{aligned}
   {d \underline{r}_{i}(t) \over d t} =& \underline{v}_{i}(t) \nonumber \\
   {d \underline{v}_{i}(t) \over d t} =& {{\underline{f}_{i}(t)+\underline{R}_{i}(t)} \over
   m_{i}} - \chi \; \underline{v}_{i}(t)~~,\end{aligned}

where :math:`\chi` is the user defined *constant* (positive, in units of
ps\ :math:`^{-1}`) specifying the thermostat friction parameter and
:math:`R(t)` is stochastic force with zero mean that satisfies the
fluctuation- dissipation theorem:

.. math::
   :label: langevin_eq

   \left< R^{\alpha}_{i}(t)~R^{\beta}_{j}(t^\prime)\right> =
   2~\chi~m_{i}~k_{B}T~\delta_{ij}~\delta_{\alpha \beta}~\delta(t-t^\prime)~~, \label{langevin}

where superscripts denote Cartesian indices, subscripts particle
indices, :math:`k_{B}` is the Boltzmann constant, :math:`T` the target
temperature and :math:`m_{i}` the particle’s mass. The Stokes-Einstein
relation for the diffusion coefficient can then be used to show that the
average value of :math:`R_{i}(t)` over a time step (in thermal
equilibrium) should be a random deviate drawn from a Gaussian
distribution of zero mean and unit variance,
:math:`\texttt{ Gauss}(0,1)`, scaled by
:math:`\sqrt{\frac{2~\chi~m_{i}~k_{B}T}{\Delta t}}`.

The effect of this algorithm is thermostat the system on a local scale.
Particles that are too “cold” are given more energy by the noise term
and particles that are too “hot” are slowed down by the friction.
Numerical instabilities, which usually arise from inaccurate calculation
of a local collision-like process, are thus efficiently kept under
control and cannot propagate.

The generation of random forces is implemented in the routine
``langevin_forces``.

An inhomogeneous variant of the Langevin thermostat
:cite:`duffy-07a` allows :math:`\chi` to vary according to
atomic velocity: it can be increased from :math:`\chi_{ep}` to
:math:`\chi_{ep} + \chi_{es}` when the atomic speed is higher than a
cut-off value. When using this thermostat as part of the two-temperature
model (TTM), :math:`\chi_{ep}` represents the required friction
parameter for electron-phonon coupling and :math:`\chi_{es}` gives the
increase due to electronic stopping. TTM simulations also require the
random deviate to be scaled by
:math:`\sqrt{\frac{2~\chi_{ep}~m_{i}~k_{B}T_{e}}{\Delta t}}`, where
:math:`T_{e}` is the local electronic temperature for the atom.

The VV implementation of the algorithm is tailored in a Langevin Impulse
(LI) manner :cite:`izaguirre-01a`:

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + \epsilon) &\leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{v}(t + {1 \over 2} \Delta t - \epsilon) &\leftarrow&
      \exp(-\chi~\Delta t) \; \underline{v}(t + \epsilon) +
      \frac{\sqrt{2~\chi~m~k_{B}T}}{m} \; \underline{Z}_{1}(\chi , \Delta t) \\
      \underline{r}(t + \Delta t) &\leftarrow& \underline{r}(t) +
      \frac{1-\exp(-\chi~\Delta t)}{\chi} \; \underline{v}(t + \epsilon) +
      \frac{\sqrt{2~\chi~m~k_{B}T}}{\chi~m} \;
      \underline{Z}_{2}(\chi , \Delta t)~~, \nonumber\end{aligned}

   where :math:`\underline{Z}_{1}(\chi , \Delta t)` and
   :math:`\underline{Z}_{2}(\chi , \Delta
   t)` are joint Gaussian random variables of zero mean, sampling from a
   *bivariate* Gaussian distribution :cite:`izaguirre-01a`:

   .. math::

      \left[
        \begin{array}{c}
          \underline{Z}_{1} \\
          \underline{Z}_{2} \\
        \end{array}
      \right]
      =
      \left[
        \begin{array}{cc}
          \sigma_{2}^{1/2} & 0 \\
          (\sigma_{1} - \sigma_{2})\sigma_{2}^{-1/2} & (\Delta t - \sigma_{1}^{2}\sigma_{2}^{-1})^{1/2} \\
        \end{array}
      \right]
      \left[
        \begin{array}{c}
          \underline{R}_{1} \\
          \underline{R}_{2} \\
        \end{array}
      \right]

   with

   .. math:: \sigma_{k} = \frac{1 - \exp(-k~\chi~\Delta t)}{k~\chi}~~,~~k=1,2

   and :math:`\underline{R}_{k}` vectors of independent standard Gaussian
   random numbers of zero mean and unit variance,
   :math:`\texttt{ Gauss}(0,1)`, - easily related to the Langevin random
   forces as defined in equation :eq:`langevin_eq`.

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t - \epsilon) +
      {\Delta t \over 2} \; {\underline{f}(t + \Delta t) \over m}

#. RATTLE_VV2  .

The algorithm is self-consistent and requires no iterations. It is worth
noting that the integration is conditional upon the Cholesky
factorisation which is impossible when
:math:`\Delta t < \sigma_{1}^{2}/\sigma_{2}` . Happily, one can show
that

.. math::

   \sigma_{1}^{2}/\sigma_{2} =
   \dfrac{2}{\chi}\tanh{\bigg(\dfrac{\chi \Delta t}{2} \bigg)},

which, for non-negative :math:`\chi, \Delta t` is always
:math:`\leq \Delta t` .

.. note::
   
   By the nature of the ensemble the centre of mass will not
   be stationary although the ensemble average warrants its proximity to
   the its original position, i.e. the COM momentum accumulation ensemble
   average will tend towards zero. By default this accumulation is removed
   and thus the correct application of stochastic dynamics the user is
   advised to use in the **no vom** option in the CONTROL file (see
   Section :ref:`control-file`). If the option is not
   applied then the dynamics will lead to peculiar thermalisation of
   different atomic species to mass- and system size-dependent
   temperatures.

The VV flavour of the Langevin thermostat is implemented in the
DL_POLY_4routines ``nvt_l0_vv``. The routines ``nvt_l1_vv`` implements
the same but also incorporate RB dynamics. The inhomogeneous Langevin
thermostat is implemented in the DL_POLY_4routine ``nvt_l2_vv`` no RB
dynamics are currently available for this form of Langevin thermostat.

Andersen Thermostat
-------------------

This thermostat assumes the idea that the system, or some subset of the
system, has an instantaneous interaction with some fictional particles
and exchanges energy. Practically, this interaction amounts to replacing
the momentum of some atoms with a new momentum drawn from the correct
Boltzmann distribution at the desired temperature. The strength of the
thermostat can be adjusted by setting the average time interval over
which the interactions occur, and by setting the magnitude of the
interaction. The collisions are best described as a random (Poisson)
process so that the probability that a collision occurs in a time step
:math:`\Delta t` is

.. math:: P_\texttt{ collision}(t) = 1 - \exp \left(-\frac{\Delta t}{\tau_{T}}\right)~~,

where :math:`\tau_{T}` is the thermostat relaxation time. The hardest
collision is to completely reset the momentum of the Poisson selected
atoms in the system, with a new one selected from the Boltzmann
distribution

.. math::

   F(\underline{v}_{i}) = \sqrt{\left(\frac{m_{i}}{2 \pi k_{B}T_\texttt{ ext}}\right)^{3}}
   \exp \left(-\frac{m_{i}~\underline{v}_{i}^{2}}{2 k_{B}T_\texttt{ ext}}\right) =
   \sqrt{\frac{k_{B}T_\texttt{ ext}}{2 m_{i}}}~\texttt{ Gauss}(0,1)~~.

where subscripts denote particle indices, :math:`k_{B}` is the Boltzmann
constant, :math:`T_\texttt{ ext}` the target temperature and
:math:`m_{i}` the particle’s mass. The thermostat can be made softer by
mixing the new momentum :math:`\underline{v}_{i}^\texttt{ new}` drawn from
:math:`F(\underline{v}_{i})` with the old momentum
:math:`\underline{v}_{i}^\texttt{ old}`

.. math::

   \underline{v}_{i} = \alpha \; \underline{v}_{i}^\texttt{ old} +
   \sqrt{1-\alpha^{2}} \; \underline{v}_{i}^\texttt{ new}~~,

where :math:`\alpha` (\ :math:`0 \le \alpha \le 1`) is the softness of the
thermostat. In practice, a uniform distribution random number,
:math:`\texttt{ uni}(i)`, is generated for each particle in the system,
which is compared to the collision probability. If
:math:`\texttt{ uni}(i) \le 1 - \exp \left(-\frac{\Delta t}{\tau_{T}}\right)`
the particle momentum is changed as described above.

The VV implementation of the Andersen algorithm is as follows:

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) &\leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) &\leftarrow& \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t)\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; \left[{\underline{f}(t + \Delta t) \over m} \right]

#. RATTLE_VV2

#. Thermostat: Note that the MD cell centre of mass momentum must not
   change!

   .. math::

      \begin{aligned}
      \texttt{ If} && \left(\texttt{ uni}(i) \le 1 - \exp \left(-\frac{\Delta t}{\tau_{T}}\right) \right)
      \; \texttt{ Then} \nonumber \\
      \underline{v}_{i}^\texttt{ new}(t + \Delta t) &\leftarrow&
      \sqrt{\frac{k_{B}T}{2 m_{i}}}~\texttt{ Gauss}(0,1) \\
      \underline{v}_{i}(t + \Delta t) &\leftarrow& \alpha \; \underline{v}_{i}(t + \Delta t) +
      \sqrt{1-\alpha^{2}} \; \underline{v}_{i}^\texttt{ new}(t + \Delta t) \nonumber \\
      \texttt{ End} && \texttt{ If}~~.\nonumber\end{aligned}

The algorithm is self-consistent and requires no iterations.

The VV flavour of the Andersen thermostat is implemented in the
DL_POLY_4routine ``nvt_a0_vv``. The routine ``nvt_a1_vv`` implements the
same but also incorporate RB dynamics.

Berendsen Thermostat
--------------------

In the Berendsen algorithm the instantaneous temperature is pushed
towards the desired temperature :math:`T_{\rm ext}` by scaling the
velocities at each step by

.. math::

   \chi (t) = \left[ 1 + {\Delta t \over \tau_{T}} \left(
   {\sigma \over E_{kin}(t)} - 1 \right) \right]^{1/2}~~,

where

.. math:: \sigma = \frac{f}{2}~k_{B}~T_{\rm ext} \label{sigma}
   :label: sigma_eq

is the target thermostat energy (depending on the external temperature
and the system total degrees of freedom, :math:`f` -
equation :eq:`freedom_eq` and :math:`\tau_{T}` a specified
time constant for temperature fluctuations (normally in the range [\ :math:`0.5`,
:math:`2``] ps).

The VV implementation of the Berendsen algorithm is straight forward. A
conventional VV1 and VV2 (thermally unconstrained) steps are carried
out. At the end of VV2 velocities are scaled by a factor of :math:`\chi`
in the following manner

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) &\leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) &\leftarrow& \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t)\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; {\underline{f}(t + \Delta t) \over m}

#. RATTLE_VV2

#. Thermostat:

   .. math::

      \begin{aligned}
      \chi (t + \Delta t) &\leftarrow& \left[ 1 + {\Delta t \over \tau_{T}}
      \left( {\sigma \over E_{kin}(t + \Delta t)} - 1 \right) \right]^{1/2} \nonumber \\
      \underline{v}(t + \Delta t) &\leftarrow& \underline{v}(t + \Delta t) \; \chi~~.\end{aligned}

.. note::
   
   The MD cell’s centre of mass momentum is removed at the
   end of the integration algorithms.

The Berendsen algorithms conserve total momentum but not energy.

The VV flavour of the Berendsen thermostat is implemented in the
DL_POLY_4routine ``nvt_b0_vv``. The routine ``nvt_b1_vv`` implements the
same but also incorporate RB dynamics.

Nosé-Hoover Thermostat
----------------------

In the Nosé-Hoover algorithm :cite:`hoover-85a` Newton’s
equations of motion are modified to read:

.. math::

   \begin{aligned}
   {d \underline{r}(t) \over d t} =& \underline{v}(t) \nonumber \\
   {d \underline{v}(t) \over d t} =& {\underline{f}(t) \over m} - \chi (t) \; \underline{v}(t)\end{aligned}

The friction coefficient, :math:`\chi`, is controlled by the first order
differential equation

.. math:: {d \chi (t) \over dt} = {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}}

where :math:`\sigma` is the target thermostat energy, equation
:eq:`sigma_eq`, and

.. math:: q_{mass} = 2~\sigma~\tau_{T}^{2}

is the thermostat mass, which depends on a specified time constant
:math:`\tau_{T}` (for temperature fluctuations normally in the range
[\ :math:`0.5`, :math:`2`] ps).

The VV implementation of the Nosé-Hoover algorithm takes place in a
symplectic manner as follows:

#. Thermostat: Note :math:`E_{kin}(t)` changes inside

   .. math::

      \begin{aligned}
      \chi (t + {1 \over 4} \Delta t) \leftarrow& \chi (t) +
      {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber \\
      \underline{v}(t) leftarrow& \underline{v}(t) \; \exp \left(
      -\chi (t + {1 \over 4} \Delta t) \; {\Delta t \over 2} \right) \\
      \chi (t + {1 \over 2} \Delta t) \leftarrow& \chi (t + {1 \over 4} \Delta t) +
      {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber \\\end{aligned}

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) \leftarrow& \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t)\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; {\underline{f}(t + \Delta t) \over m}

#. RATTLE_VV2

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` changes inside

   .. math::

      \begin{aligned}
      \chi (t + {3 \over 4} \Delta t) \leftarrow& \chi (t + {1 \over 2} \Delta t) +
      {\Delta t \over 4} \; {{2 E_{kin}(t + \Delta t) - 2 \sigma} \over q_{mass}} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + \Delta t) \; \exp \left(
      -\chi (t + {3 \over 4} \Delta t) \; {\Delta t \over 2} \right) \\
      \chi (t + \Delta t) \leftarrow& \chi (t + {3 \over 4} \Delta t) +
      {\Delta t \over 4} \; {{2 E_{kin}(t + \Delta t) - 2 \sigma} \over q_{mass}}~~. \nonumber\end{aligned}

The algorithm is self-consistent and requires no iterations.

The conserved quantity is derived from the extended Hamiltonian for the
system which, to within a constant, is the Helmholtz free energy:

.. math::

   {\cal H}_{\rm NVT} = {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   f~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~,

where :math:`f` is the system’s degrees of freedom - equation
:eq:`freedom_eq`.

The VV flavour of the Nosé-Hoover thermostat is implemented in the
DL_POLY_4routine ``nvt_h0_vv``. The routine ``nvt_h1_vv`` implements the
same but also incorporate RB dynamics.

Gentle Stochastic Thermostat
----------------------------

The Gentle Stochastic Thermostat
:cite:`leimkuhler-09a,samoletov-07a` is an extension of the
Nosé-Hoover algorithm :cite:`hoover-85a`

.. math::

   \begin{aligned}
   {d \underline{r}(t) \over d t} =& \underline{v}(t) \nonumber \\
   {d \underline{v}(t) \over d t} =& {\underline{f}(t) \over m} - \chi (t) \; \underline{v}(t)\end{aligned}

in which the thermostat friction, :math:`\chi`, has its own Brownian
dynamics:

.. math::
   :label: Ornstein-Uhlenbeck_eq

   {d \chi (t) \over dt} = \frac{2 E_{kin}(t) - 2 \sigma}{q_{mass}} - \gamma~\chi(t) +
   \frac{\sqrt{2~\gamma~k_{B}~T_{\rm ext}~q_{mass}}}{q_{mass}}~{d \omega(t) \over d t}~~, \label{Ornstein-Uhlenbeck}

governed by the Langevin friction :math:`\gamma` (positive, in units of
ps\ :math:`^{-1}`), where :math:`\omega(t)` is the standard Brownian
motion (Wiener process - ``Gauss``\ (0,1)), :math:`\sigma` is the target
thermostat energy, as in equation :eq:`sigma_eq`.

.. math:: q_{mass} = 2~\sigma~\tau_{T}^{2}

is the thermostat mass, which depends on a specified time constant
:math:`\tau_{T}` (for temperature fluctuations normally in the range
[\ :math:`0.5`, :math:`2`] ps).

It is worth noting that
equation :eq:`Ornstein-Uhlenbeck_eq` similar to
the Ornstein-Uhlenbeck equation:

.. math:: {d \chi \over dt} = -\frac{\alpha\sigma^{2}}{2}\chi + \sigma {d \omega \over dt}~~,

which for a given realization of the Wiener process :math:`\omega(t)`
has an exact solution:

.. math::

   \chi_{n+1} = e^{-\epsilon t} \left( \chi_{n} + \sigma
   \sqrt{\frac{e^{2~\epsilon t} - 1}{2 \epsilon}} \Delta \omega \right)~~,

where :math:`\epsilon = \alpha\sigma^{2}/2` and
:math:`\Delta \omega \sim \mathcal{N}(0,1)`. The VV implementation of
the Gentle Stochastic Thermostat algorithm takes place in a symplectic
manner as follows:

#. Thermostat: Note :math:`E_{kin}(t)` changes inside and
   :math:`R_{g_{1,2}}(t)`, drawn from :math:`\texttt{ Gauss}(0,1)`, are
   independent

   .. math::

      \begin{aligned}
      \chi (t + {1 \over 4} \Delta t) \leftarrow& \chi (t) ~ \exp \left[-\gamma~{\Delta t \over 4} \right] +
      \sqrt{\frac{k_{B}~T_{\rm ext}}{q_{mass}}\left(1-\exp^{2} \left[-\gamma~{\Delta t \over 4} \right]\right)}~R_{g_{1}}(t) \phantom{xxx} \nonumber \\
      &\phantom{xxxxxxxxxxxxxxxx} + {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber \\
      \underline{v}(t) \leftarrow& \underline{v}(t) \; \exp \left(
      -\chi (t + {1 \over 4} \Delta t) \; {\Delta t \over 2} \right) \\
      \chi (t + {1 \over 2} \Delta t) \leftarrow& \chi (t + {1 \over 4} \Delta t) ~ \exp \left[-\gamma~{\Delta t \over 4} \right] \nonumber \\
      &~~~+\sqrt{\frac{k_{B}~T_{\rm ext}}{q_{mass}}\left(1-\exp^{2} \left[-\gamma~{\Delta t \over 4} \right]\right)}~R_{g_{2}}(t+{1 \over 4} \Delta t)
      + {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber\end{aligned}

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) \leftarrow& \underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t)\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \underline{v}(t + \Delta t) \leftarrow \underline{v}(t + {1 \over 2} \Delta t) +
      {\Delta t \over 2} \; {\underline{f}(t + \Delta t) \over m}

#. RATTLE_VV2

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` changes inside and
   :math:`R_{g_{3,4}}(t)`, drawn from :math:`\texttt{ Gauss}(0,1)`, are
   independent

   .. math::

      \begin{aligned}
      \chi (t + {3 \over 4} \Delta t) \leftarrow& \chi (t + {1 \over 2} \Delta t) ~ \exp \left[-\gamma~{\Delta t \over 4} \right] +
      \sqrt{\frac{k_{B}~T_{\rm ext}}{q_{mass}}\left(1-\exp^{2} \left[-\gamma~{\Delta t \over 4} \right]\right)}~R_{g_{3}}(t + {2 \over 4} \Delta t) \phantom{xxx} \nonumber \\
      &\phantom{xxxxxxxxxxxxxxxx} + {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + \Delta t) \; \exp \left(
      -\chi (t + {3 \over 4} \Delta t) \; {\Delta t \over 2} \right) \\
      \chi (t + \Delta t) \leftarrow& \chi (t + {3 \over 4} \Delta t) ~ \exp \left[-\gamma~{\Delta t \over 4} \right] +
      \sqrt{\frac{k_{B}~T_{\rm ext}}{q_{mass}}\left(1-\exp^{2} \left[-\gamma~{\Delta t \over 4} \right]\right)}~R_{g_{4}}(t + {3 \over 4} \Delta t) \phantom{xxx} \nonumber \\
      &\phantom{xxxxxxxxxxxxxxxx} + {\Delta t \over 4} \; {{2 E_{kin}(t) - 2 \sigma} \over q_{mass}} \nonumber\end{aligned}

The algorithm is self-consistent and requires no iterations.

The conserved quantity is derived from the extended Hamiltonian for the
system which, to within a constant, is the Helmholtz free energy:

.. math::

   {\cal H}_{\rm NVT} = {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   f~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~,

where :math:`f` is the system’s degrees of freedom - equation
:eq:`freedom_eq`.

The VV flavour of the Gentle Stochastic Thermostat is implemented in the
DL_POLY_4routine ``nvt_g0_vv``. The routine ``nvt_g1_vv`` implements the
same but also incorporate RB dynamics.

.. _dpd:

Dissipative Particle Dynamics Thermostat
----------------------------------------

An elegant way to integrate the DPD equations of motions, as shown in
:ref:`Appendix A <DPD-all>`, is introduced by Shardlow
:cite:`shardlow-03a`. By applying ideas commonly used in
solving differential equations to the case of integrating the equations
of motion in DPD, the integration process is factorised by splitting the
conservative forces calculation from that of the dissipative and random
terms. In this way the conservative part can be solved using traditional
molecular dynamics methods, while the fluctuation-dissipation part is
solved separately as a stochastic differential (Langevin) equation.
There are two Shardlow integrators, called S1 ( **dpds1**) and S2 (
**dpds2**), based on splitting the equations of motion up to first and
second order, respectively, using Suzuki-Trotter(Strang) expansion of
the Liouville evolution operator and thus warranting the integrators’
symplectic.

To describe the two integrators we define the algorithmic sequence

.. math::

   {\rm S}(\Delta t)~~=~~\left\{ \begin{array} {l}
   \textnormal{~For all pairs of particles for which}~r_{ij} < r_{c} : \\
   \begin{array} {l}
   (1):~~~\underline{v}_{i}~\leftarrow~\underline{v}_{i} - \frac{1}{2 m_{i}} \left\{ \gamma_{ij} w^{2}_{(r_{ij})}
   (\underline{v}_{ij} \cdot \underline{e}_{ij}) \underline{e}_{ij} \Delta t +
   \sqrt{2 \gamma_{ij} k_{B} T} w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} \right\} \\
   (2):~~~\underline{v}_{j}~\leftarrow~\underline{v}_{j} + \frac{1}{2 m_{i}} \left\{ \gamma_{ij} w^{2}_{(r_{ij})}
   (\underline{v}_{ij} \cdot \underline{e}_{ij}) \underline{e}_{ij} \Delta t -
   \sqrt{2 \gamma_{ij} k_{B} T} w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} \right\} \\
   (3):~~~\underline{v}_{i}~\leftarrow~\underline{v}_{i} + \frac{1}{2 m_{i}} \left\{ \sqrt{2 \gamma_{ij} k_{B} T}
   w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} -
   \frac{\gamma_{ij} w^{2}(r_{ij}) \Delta t}{1 + \gamma_{ij} w^{2}(r_{ij}) \Delta t} \right. \times \\
   \phantom{xxxxxxxxxxxxxxxxxxxxxxxxx} \left. \left[(\underline{v}_{ij} \cdot \underline{e}_{ij}) \underline{e}_{ij} +
   \sqrt{2 \gamma_{ij} k_{B} T} w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} \right] \right\} \\
   (4):~~~\underline{v}_{j}~\leftarrow~\underline{v}_{j} - \frac{1}{2 m_{i}} \left\{ \sqrt{2 \gamma_{ij} k_{B} T}
   w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} +
   \frac{\gamma_{ij} w^{2}(r_{ij}) \Delta t}{1 + \gamma_{ij} w^{2}(r_{ij}) \Delta t}  \right. \times \\
   \phantom{xxxxxxxxxxxxxxxxxxxxxxxxx} \left. \left[(\underline{v}_{ij} \cdot \underline{e}_{ij}) \underline{e}_{ij} +
   \sqrt{2 \gamma_{ij} k_{B} T} w_{(r_{ij})} \zeta_{ij} \underline{e}_{ij} \sqrt{\Delta t} \right] \right\}
   \end{array} \end{array} \right.

as the Shardlow operator, where :math:`\zeta_{ij}` is the random number
with zero mean and unit variance, unique for every unique pair
:math:`\left\{ij\right\}` in the system, and
:math:`w_{(r_{ij})} = w^{C}(r_{ij})` is the DPD conservative force
switching function in equation :eq:`DPDS_eq`,
:math:`\gamma_{ij}` is the drag coefficient between the types of
particles :math:`i` and :math:`j`,
:math:`\underline{v}_{ij} = \underline{v}_{j}-\underline{v}_{i}` is the inter-particle
relative velocity and :math:`\underline{e}_{ij} = \underline{r}_{ij}/r_{ij}` is the
inter-particle unit vector. :math:`k_{B}` and :math:`T` are the
Boltzmann constant and the system target temperature.

If we define the velocity Verlet micro-canonical (NVE) ensemble sequence
as :math:`{\rm NVE}(t + \Delta t)` then Shardlow’s first
(:math:`{\cal S}1`) and second (:math:`{\cal S}2`) order splittings can
be written algorithmically as the following sequential operators
applications:

.. math::

   \begin{aligned}
   {\cal S}1 ::& {\rm S(~\Delta t~~)} \rightarrow {\rm NVE}(t + \Delta t) \nonumber \\
   {\cal S}2 ::& {\rm S}(\Delta t / 2) \rightarrow {\rm NVE}(t + \Delta t) \rightarrow {\rm S}(\Delta t / 2)~~.\end{aligned}

The application of these DPD thermostats are implemented in the
DL_POLY_4routine ``dpd_thermostat`` and which only applies as a
perturbation around the NVE integrator incorporating both particle and
RB dynamics.
