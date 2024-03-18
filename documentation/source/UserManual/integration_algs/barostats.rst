Barostats
=========

The size and shape of the simulation cell may be dynamically adjusted by
coupling the system to a barostat in order to obtain a desired average
pressure (:math:`P_{\rm ext}`) and/or isotropic :index:`stress tensor`
(:math:`\underline{\underline{\mathbf{\sigma}}}`). DL_POLY_4has four such algorithms: the Langevin
type barostat :cite:`quigley-04a`, the Berendsen barostat
:cite:`berendsen-84a`, the Nosé-Hoover type barostat
:cite:`hoover-85a` and the Martyna-Tuckerman-Klein (MTK)
barsotat :cite:`martyna-96a`. Only the Berendsen barostat
does not have defined conserved quantity.

.. note::
   
   The MD cell’s centre of mass momentum is removed at the
   end of the integration algorithms with barostats.

Instantaneous pressure and stress
---------------------------------

The instantaneous pressure in a system,

.. math::
   :label: prs_inst_eq

   {\cal P}(t) = \frac{\left[ 2 E_{kin}(t) - {\cal W}_{\rm atomic}(t) -
   {\cal W}_{\rm constrain}(t - \Delta t) -
   {\cal W}_{\rm PMF}(t - \Delta t) \right]} {3 V(t)}~~, \label{prs_inst}

is a function of the system volume, kinetic energy and virial,
:math:`{\cal W}`. 

.. note::
   
   when bond constraints or/and PMF
   constraints are present in the system :math:`{\cal P}` will not converge
   to the exact value of :math:`P_{\rm ext}` during equilibration in NPT
   and N\ :math:`\sigma`\ T simulations. This is due to iterative nature of
   the constrained motion, in which the virials
   :math:`{\cal W}_{\rm constrain}` and :math:`{\cal W}_{\rm PMF}` are
   calculated retrospectively to the forcefield virial
   :math:`{\cal W}_{\rm atomic}`.

The instantaneous stress tensor in a system,

.. math::
   :label: str_inst_eq

   \underline{\underline{\mathbf{\sigma}}}(t) = \underline{\underline{\mathbf{\sigma}}}_{kin}(t) + \underline{\underline{\mathbf{\sigma}}}_{\rm atomic}(t) +
   \underline{\underline{\mathbf{\sigma}}}_{\rm constrain}(t - \Delta t) + \underline{\underline{\mathbf{\sigma}}}_{\rm PMF}(t - \Delta t)~~, \label{str_inst}

is a sum of the forcefield, :math:`\underline{\underline{\mathbf{\sigma}}}_{\rm atomic}`,
constrain, :math:`\underline{\underline{\mathbf{\sigma}}}_{\rm constrains}`, and PMF,
:math:`\underline{\underline{\mathbf{\sigma}}}_{\rm PMF}`, stresses. 

.. note::
   
   When bond
   constraints or/and PMF constraints are present in the system, the
   quantity :math:`\frac{\texttt{ Tr}[\underline{\underline{\mathbf{\sigma}}}]}{3 V}` will not
   converge to the exact value of :math:`P_{\rm ext}` during equilibration
   in NPT and N\ :math:`\sigma`\ T simulations. This is due to iterative
   nature of the constrained motion in which the constraint and PMF
   stresses are calculated retrospectively to the forcefield stress.

Langevin Barostat
-----------------

DL_POLY_4implements a Langevin barostat :cite:`quigley-04a`
for isotropic and anisotropic cell fluctuations.

Cell size variations
~~~~~~~~~~~~~~~~~~~~

For isotropic fluctuations the equations of motion are:

.. math::

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \eta (t) \; \underline{r}(t) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t) + \underline{R}(t)}{m} - \left[ \chi +
   \left(1+\frac{3}{f}\right) \eta (t) \right] \underline{v}(t) \nonumber \\
   \frac{d}{dt}\eta (t) =& 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} +
   3 \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} - \chi_{p}~\eta (t) + \frac{R_{p}}{p_{mass}} \\
   p_{mass} =& \frac{(f+3)~k_{B}~T_{\rm ext}}{(2 \pi~\chi_{p})^{2}} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \eta (t) \; \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& [3 \eta (t)]~V(t)~~, \nonumber\end{aligned}

where :math:`\chi` and :math:`\chi_{p}` are the user defined *constants*
(positive, in units of ps\ :math:`^{-1}`), specifying the thermostat and
barostat friction parameters, :math:`R(t)` is the Langevin stochastic
force (see equation :eq:`langevin_eq`), :math:`{\cal P}`
the instantaneous pressure (equation :eq:`prs_inst_eq`) and
:math:`R_{p}` is the stochastic (Langevin) pressure variable

.. math:: \left< R_{p}(t)~R_{p}(t^\prime)\right> = 2~\chi_{p}~p_{mass}~k_{B}T~\delta(t-t^\prime)~~,

which is drawn from Gaussian distribution of zero mean and unit
variance, :math:`\texttt{ Gauss}(0,1)`, scaled by //
:math:`\sqrt{\frac{2~\chi_{p}~p_{mass}~k_{B}T}{\Delta t}}`.
:math:`k_{B}` is the Boltzmann constant, :math:`T` the target
temperature and :math:`p_{mass}` the barostat mass. is the cell matrix
whose columns are the three cell vectors
:math:`\underline{a}, \underline{b}, \underline{c}`.

The conserved quantity these generate is:

.. math:: {\cal H}_{\rm NPT} = {\cal H}_{\rm NVE} + {p_{mass}~\eta (t)^{2} \over 2} + P_{\rm ext} V(t)~~.

The VV implementation of the Langevin algorithm only requires iterations
if bond or PMF constraints are present (:math:`4` until satisfactory
convergence of the constraint forces is achieved). These are with
respect to the pressure (i.e. :math:`\eta (t)`) in the first part,
VV1+RATTLE_VV1. The second part is conventional, VV2+RATTLE_VV2, as at
the end the velocities are scaled by a factor of :math:`\chi`.

#. Thermostat: Note :math:`E_{kin}(t)` changes inside

   .. math:: \underline{v}(t) \leftarrow \exp \left( -\chi \; {\Delta t \over 4} \right) \; \underline{v}(t)

#. Barostat: Note :math:`E_{kin}(t)` and :math:`{\cal P}(t)` have
   changed and change inside

   .. math::

      \begin{aligned}
      \eta (t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right) \;
      \eta (t)\nonumber \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \eta (t) + {\Delta t \over 4} \;
      \left[ 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} + \right. \nonumber \\
      & ~~~~~~~~~~~~~~~~~~~~~~~~\left. 3 \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} + \frac{R_{p}(t)}{p_{mass}} \right] \nonumber \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + {1 \over 4} \Delta t) \nonumber \\
      \underline{v}(t) \leftarrow& \exp \left[ -\left( 1 + \frac{3}{f} \right)
      \eta (t + {1 \over 4}\Delta t) \; {\Delta t \over 2} \right] \; \underline{v}(t) \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + {1 \over 4} \Delta t) \nonumber \\
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \eta (t + {1 \over 4} \Delta t) + {\Delta t \over 4} \;
      \left[ 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} + \right. \nonumber \\
      & ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\left. 3 \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} + \frac{R_{p}(t)}{p_{mass}} \right] \nonumber \\
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + {1 \over 2} \Delta t) \nonumber\end{aligned}

#. Thermostat: Note :math:`E_{kin}(t)` has changed and changes inside

   .. math:: \underline{v}(t) \leftarrow \exp \left( -\chi \; {\Delta t \over 4} \right) \; \underline{v}(t)

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; \frac{\underline{f}(t)+\underline{R}(t)}{m} \nonumber \\
      \underline{\underline{\mathbf{H}}}(t + \Delta t) \leftarrow& \exp \left[
      \eta (t + {1 \over 2} \Delta t) \; \Delta t \right] \; \underline{\underline{\mathbf{H}}}(t) \nonumber \\
      V(t + \Delta t) \leftarrow& \exp \left[3 \eta (t + {1 \over 2} \Delta t) \;
      \Delta t \right] \; V(t) \\
      \underline{r}(t + \Delta t) \leftarrow& \exp \left[ \eta (t + {1 \over 2} \Delta t) \; \Delta t \right] \;
      \underline{r}(t) + \Delta t \; \underline{v}(t + {1 \over 2} \Delta t) \nonumber\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math::

      \begin{aligned}
      \underline{f}(t + \Delta t) \leftarrow& \underline{f}(t) \nonumber \\
      \underline{R}(t + \Delta t) \leftarrow& \underline{R}(t) \\
      R_{p} (t + \Delta t) \leftarrow& R_{p} (t) \nonumber\end{aligned}

#. VV2:

   .. math::

      \begin{aligned}
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + {\Delta t \over 2}) +
      {\Delta t \over 2} \; \frac{\underline{f}(t)+\underline{R}(t)}{m}\end{aligned}

#. RATTLE_VV2

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` has changed and
   changes inside

   .. math:: \underline{v}(t + \Delta t) \leftarrow \exp \left( -\chi \; {\Delta t \over 4} \right) \; \underline{v}(t + \Delta t)

#. Barostat: Note :math:`E_{kin}(t + \Delta t)` and
   :math:`{\cal P}(t + \Delta t)` have changed and change inside

   .. math::

      \begin{aligned}
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right) \;
      \eta (t + {1 \over 2} \Delta t) \nonumber \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \eta (t + {1 \over 2} \Delta t) + {\Delta t \over 4} \;
      \left[ 3 V(t + \Delta t) \frac{{\cal P}(t + \Delta t) - P_{\rm ext}}{p_{mass}} + \right. \nonumber \\
      & ~~~~~~~~~~~~~~~~~~~~~~~~~~~\left. 3 \frac{2 E_{kin}(t + \Delta t)}{f} \frac{1}{p_{mass}} + \frac{R_{p}(t)}{p_{mass}} \right] \nonumber \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + {3 \over 4} \Delta t) \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \exp \left[ -\left( 1 + \frac{3}{f} \right)
      \eta (t + {3 \over 4} \Delta t) \; {\Delta t \over 2} \right] \; \underline{v}(t + \Delta t) \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + {3 \over 4} \Delta t) \nonumber \\
      \eta (t + \Delta t) \leftarrow& \eta (t + {3 \over 4} \Delta t) + {\Delta t \over 4} \;
      \left[ 3 V(t + \Delta t) \frac{{\cal P}(t + \Delta t) - P_{\rm ext}}{p_{mass}} + \right. \nonumber \\
      & ~~~~~~~~~~~~~~~~~~~~~~~~~~~\left. 3 \frac{2 E_{kin}(t + \Delta t)}{f} \frac{1}{p_{mass}} + \frac{R_{p}(t)}{p_{mass}} \right] \nonumber \\
      \eta (t + \Delta t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right)  \;
      \eta (t + \Delta t) \nonumber\end{aligned}

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` has changed and
   changes inside

   .. math:: \underline{v}(t + \Delta t) \leftarrow \exp \left( -\chi \; {\Delta t \over 4} \right) \; \underline{v}(t + \Delta t)~~,

The VV flavour of the langevin barostat (and Nosé-Hoover thermostat) is
implemented in the DL_POLY_4routine ``npt_l0_vv``. The routine
``npt_l1_vv`` implements the same but also incorporate RB dynamics.

Cell size and shape variations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The isotropic algorithms may be extended to allowing the cell shape to
vary by defining :math:`\eta` as a tensor, :math:`\underline{\underline{\mathbf{\eta}}}` and
extending the Langevin pressure variable :math:`R_{p}` to a stochastic
(Langevin) tensor :math:`\underline{\underline{\mathbf{R_{p}}}}`:

.. math:: \left< R_{p,i}(t)~R_{p,j}(t^\prime)\right> = 2~\chi_{p}~p_{mass}~k_{B}T~\delta_{ij}~\delta(t-t^\prime)~~,

which is drawn from Gaussian distribution of zero mean and unit
variance, :math:`\texttt{ Gauss}(0,1)`, scaled by
:math:`\sqrt{\frac{2~\chi_{p}~p_{mass}~k_{B}T}{\Delta t}}`.
:math:`k_{B}` is the Boltzmann constant, :math:`T` the target
temperature and :math:`p_{mass}` the barostat mass. **Note** that
:math:`\underline{\underline{\mathbf{R_{p}}}}` has to be symmetric and only 6 independent
components must be generated each timestep.

The equations of motion are written in the same fashion as is in the
isotropic algorithm with slight modifications (as now the equations with
:math:`\eta` are extended to matrix forms)

.. math::

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \underline{\underline{\mathbf{\eta (t)}}} \cdot \underline{r}(t) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t) + \underline{R}(t)}{m} - \left[ \chi~\underline{\underline{\mathbf{1}}} +
   \underline{\underline{\mathbf{\eta}}}(t) + \frac{\texttt{ Tr}\left[\underline{\underline{\mathbf{\eta}}}(t)\right]}{f}~\underline{\underline{\mathbf{1}}} \right] \cdot \underline{v}(t) \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{\eta}}}(t) =& \frac{\underline{\underline{\mathbf{\sigma}}}(t) -
   P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}} + \frac{2 E_{kin}(t)}{f} \frac{\underline{\underline{\mathbf{1}}}}{p_{mass}} -
   \chi_{p} \underline{\underline{\mathbf{\eta}}}(t)  + \frac{\underline{\underline{\mathbf{R_{p}}}}}{p_{mass}} \\
   p_{mass} =& \frac{(f+3)}{3}~\frac{k_{B}~T_{\rm ext}}{(2 \pi~\chi_{P})^{2}} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \underline{\underline{\mathbf{\eta}}}(t) \cdot \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& \texttt{ Tr} [\underline{\underline{\mathbf{\eta}}}(t)]~V(t)~~. \nonumber\end{aligned}

where :math:`\underline{\underline{\mathbf{\sigma}}}` is the stress tensor
(equation :eq:`str_inst_eq`) and :math:`\underline{\underline{\mathbf{1}}}` is the
identity matrix.

The conserved quantity these generate is:

.. math:: {\cal H}_{\rm N\underline{\underline{\mathbf{\sigma}}}T} = {\cal H}_{\rm NVE} + {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} + P_{\rm ext} V(t)~~.

the VV algorithmic equations are, therefore, written in the same fashion
as in the isotropic case with slight modifications. For the VV couched
algorithm these are of the following sort

.. math::

   \begin{aligned}
   \underline{\underline{\mathbf{\eta}}} (t) \leftarrow& \exp \left( -\chi_{p} \; {\Delta t \over 8} \right) \;
   \underline{\underline{\mathbf{\eta}}} (t)\nonumber \\
   \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4} \Delta t) \leftarrow& \underline{\underline{\mathbf{\eta}}}(t) + \nonumber \\
   & {\Delta t \over 4} \; \left[ \frac{\underline{\underline{\mathbf{\sigma}}}(t) - P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f} \frac{\underline{\underline{\mathbf{1}}}}{p_{mass}} + \frac{\underline{\underline{\mathbf{R_{p}}}}(t)}{p_{mass}} \right] \\
   \underline{v}(t) \leftarrow& \exp \left[ -\left( \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4}\Delta t) +
   \frac{1}{f} \texttt{ Tr}\left[\underline{\underline{\mathbf{\eta}}}(t + {1 \over 4}\Delta t)\right] \right) \;
   {\Delta t \over 2} \right] \cdot \underline{v}(t) \nonumber \\
   \underline{r}(t + \Delta t) \leftarrow& \exp \left[ \underline{\underline{\mathbf{\eta}}} (t + {1 \over 2} \Delta t) \; \Delta t \right] \cdot
   \underline{r}(t) + \Delta t \; \underline{v}(t + {1 \over 2} \Delta t) \nonumber\end{aligned}

This ensemble is optionally extending to constant normal pressure and
constant surface area, NP\ :math:`_{n}`\ AT
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion to:

.. math::

   \frac{d}{dt} \eta_{\alpha\beta}(t) = \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} + \frac{2 E_{kin}(t)}{f~p_{mass}} -
   \chi_{p} \eta_{zz}(t) + \frac{R_{p,zz}(t)}{p_{mass}} & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha,\beta) \ne z~~.
   \end{array} \right.

Similarly, this ensemble is optionally extending to constant normal
pressure and constant surface tension, NP\ :math:`_{n}\gamma`\ T
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion to:

.. math::

   \frac{d}{dt} \eta_{\alpha\beta}(t) = \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{\alpha\alpha}(t) - \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f~p_{mass}} - \chi_{p} \eta_{\alpha\alpha}(t) +
   \frac{R_{p,\alpha\alpha}(t)}{p_{mass}} & (\alpha = \beta) = x,y \\
   & \\
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} - \chi_{p} \eta_{zz}(t) +
   \frac{R_{p,zz}(t)}{p_{mass}} & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha \ne \beta) = x,y,z
   \end{array} \right. ,

where :math:`\gamma_{\rm ext}` is the user defined external surface
tension and :math:`h_{z}(t) = V(t) / A_{xy}(t)` is the instantaneous
hight of the MD box (or MD box volume over area). The instnatneous
surface tension is defined as

.. math:: \gamma_{\alpha}(t)=-h_{z}(t)\left[ \sigma_{\alpha\alpha}(t) - P_{\rm ext} \right]~~.\label{gamma}
   :label: gamma_eq

The case :math:`\gamma_{\rm ext}=0` generates the NPT anisotropic
ensemble for the orthorhombic cell (``imcon``\ :math:`=2` in CONFIG, see
:ref:`Appendix B<boundary-conditions>`). This can
be considered as an "orthorhombic" constraint on the
N\ :math:`\sigma`\ T ensemble. The constraint can be strengthened
further, to a "semi-orthorhombic" one, by imposing that the MD cell
change isotropically in the :math:`(x,y)` plane which leads to the
following modification in the N\ :math:`P_{n}\gamma`\ T set of equatons

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\alpha}(t) = \frac{\left[\sigma_{xx}(t)+\sigma_{yy}(t)\right]/2 -
   \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} +
   \frac{2~E_{kin}(t)}{f~p_{mass}} - \\
   \phantom{xxxxxxxxxxxxxx}- \chi_{p} \eta_{\alpha\alpha}(t) +
   \frac{R_{p,xx}(t)+R_{p,yy}(t)}{2~p_{mass}}~~:~~(\alpha = \beta) = x,y~~.\nonumber\end{aligned}

The VV flavour of the non-isotropic Langevin barostat (and Nosé-Hoover
thermostat) is implemented in the DL_POLY_4routine ``nst_l0_vv``. The
routine ``nst_l1_vv`` implements the same but also incorporate RB
dynamics.


.. index:: single: barostat;Berendsen

Berendsen Barostat
------------------

With the Berendsen barostat the system is made to obey the equation of
motion at the beginning of each step

.. math:: {d{\cal P}(t) \over dt} = {{P_{\rm ext} - {\cal P}(t)} \over \tau_{P}}~~,

where :math:`{\cal P}` is the instantaneous pressure
(equation :eq:`prs_inst_eq`) and :math:`\tau_{P}` is the
barostat relaxation time constant.

Cell size variations
~~~~~~~~~~~~~~~~~~~~

In the isotropic implementation, at each step the MD cell volume is
scaled by a factor :math:`\eta`, and the coordinates and cell vectors by
:math:`\eta^{1/3}`,

.. math::
   :label: berbar_eq

   \eta (t) = 1 - {\beta \Delta t \over \tau_{P}} \; (P_{\rm ext} -
   {\cal P}(t)) \label{berbar}

where :math:`\beta` is the isothermal compressibility of the system. In
practice :math:`\beta` is a specified constant which DL_POLY_4takes to
be the isothermal compressibility of liquid water. The exact value is
not critical to the algorithm as it relies on the ratio
:math:`\tau_{P}/\beta`. :math:`\tau_{P}` is a specified time constant
for pressure fluctuations, supplied by the user.

It is worth noting that the barostat and the thermostat are independent
and fully separable.

The VV implementation of the Berendsen algorithm only requires
iterations if bond or PMF constraints are present (:math:`13` until
satisfactory convergence of the constraint forces is achieved). These
are with respect to the pressure (i.e. :math:`\eta (t)`) in the first
part, VV1+RATTLE_VV1. The second part is conventional, VV2+RATTLE_VV2,
as at the end the velocities are scaled by a factor of :math:`\chi`.

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{r}(t + \Delta t) \leftarrow& \eta (t)^{1/3}~\underline{r}(t) + \Delta t \;
      \underline{v}(t + {1 \over 2} \Delta t) \\
      \underline{\underline{\mathbf{H}}}(t + \Delta t) \leftarrow&  \eta (t)^{1/3}~\underline{\underline{\mathbf{H}}}(t) \nonumber \\
      V(t + \Delta t) \leftarrow& \eta (t)~V(t) \nonumber\end{aligned}

#. RATTLE_VV1

#. Barostat:

   .. math::

      \eta (t) = 1 - {\beta \Delta t \over \tau_{P}} \; (P_{\rm ext} -
      {\cal P}(t))

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
      \chi (t + \Delta t) \leftarrow& \left[ 1 + {\Delta t \over \tau_{T}}
      \left( {\sigma \over E_{kin}(t + \Delta t)} - 1 \right) \right]^{1/2} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + \Delta t) \; \chi~~.\end{aligned}

where is the cell matrix whose columns are the three cell vectors
:math:`\underline{a}, \underline{b}, \underline{c}`.

The Berendsen algorithms conserve total momentum but not energy.

The VV flavour of the Berendsen barostat (and thermostat) is implemented
in the DL_POLY_4routine ``npt_b0_vv``. The routines ``npt_b1_vv``
implements the same but also incorporate RB dynamics.

Cell size and shape variations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The extension of the isotropic algorithm to anisotropic cell variations
is straightforward. A tensor is defined as

.. math::

   \underline{\underline{\mathbf{\eta}}}(t) = \underline{\underline{\mathbf{1}}} - {\beta \Delta t \over \tau_{P}}
   (P_{\rm ext}~\underline{\underline{\mathbf{1}}} - \underline{\underline{\mathbf{\sigma}}}(t) / V(t))~~,

where where :math:`\underline{\underline{\mathbf{\sigma}}}` is the stress tensor
(equation :eq:`str_inst_eq`) and :math:`\underline{\underline{\mathbf{1}}}` is the
identity matrix. Then new cell vectors and volume are given by

.. math::

   \begin{aligned}
   \underline{\underline{\mathbf{H}}}(t + \Delta t) \leftarrow& \underline{\underline{\mathbf{\eta}}}(t) \cdot \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   V(t + \Delta t) \leftarrow& \texttt{ Tr} [\underline{\underline{\mathbf{\eta}}}(t)]~V(t)~~.\end{aligned}

and the velocity updates as

.. math::

   \begin{aligned}
   \texttt{ VV1:}~~\underline{r}(t + \Delta t) \leftarrow& \underline{\underline{\mathbf{\eta}}}(t) \cdot \underline{r}(t) + \Delta t \;
   \underline{v}(t + {1 \over 2} \Delta t) \nonumber\end{aligned}

This ensemble is optionally extending to constant normal pressure and
constant surface area, NP\ :math:`_{n}`\ AT
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion to:

.. math::

   \eta_{\alpha\delta}(t) = \left\{ \begin{array} {l@{\quad:\quad}l}
   1 - {\beta \Delta t \over \tau_{P}} \left[ P_{\rm ext} - \sigma_{zz}(t) / V(t) \right]
   & (\alpha = \delta) = z \\
   1 & (\alpha = \delta) = x,y \\
   0 & (\alpha \ne \delta)~~.
   \end{array} \right.

Similarly, this ensemble is optionally extending to constant normal
pressure and constant surface tension, NP\ :math:`_{n}\gamma`\ T
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion to:

.. math::

   \eta_{\alpha\delta}(t) = \left\{ \begin{array} {l@{\quad:\quad}l}
   1 - {\beta \Delta t \over \tau_{P}} \left[ P_{\rm ext} - \gamma_{\rm ext}~V(t) / h_{z}(t) -
   \sigma_{\alpha\alpha}(t) / V(t) \right] & (\alpha = \delta) = x,y \\
   & \\
   1 - {\beta \Delta t \over \tau_{P}} \left[ P_{\rm ext} -
   \sigma_{zz}(t) / V(t) \right] & (\alpha = \delta) = z \\
   0 & (\alpha \ne \delta)~~,
   \end{array} \right.

where :math:`\gamma_{\rm ext}` is the user defined external surface
tension and :math:`h_{z}(t) = V(t) / A_{xy}(t)` is the instantaneous
hight of the MD box (or MD box volume over area). One defines the
instantaneous surface tension as given in
equation :eq:`gamma_eq`. The case :math:`\gamma_{\rm ext}=0`
generates the NPT anisotropic ensemble for the orthorhombic cell
(imcon=2 in CONFIG, see
:ref:`Appendix B<boundary-conditions>`). This can
be considered as an "orthorhombic" constraint on the
N\ :math:`\sigma`\ T ensemble. The constraint can be strengthened
further, to a "semi-orthorhombic" one, by imposing that the MD cell
change isotropically in the :math:`(x,y)` plane which leads to the
following change in the equations above

.. math::

   \eta_{\alpha\alpha}(t) = 1 - {\beta \Delta t \over \tau_{P}}
   \left[ P_{\rm ext} - \gamma_{\rm ext}~ \frac{V(t)}{h_{z}(t)} -
   \frac{\sigma_{xx}(t)+\sigma_{yy}(t)}{2~V(t)} \right]~~:~~(\alpha = \delta) = x,y~~.

The VV flavour of the non-isotropic Berendsen barostat (and thermostat)
is implemented in the DL_POLY_4routine ``nst_b0_vv``. The routine
``nst_b1_vv`` implements the same but also incorporate RB dynamics.


.. index:: single: barostat;Nosé-Hoover

Nosé-Hoover Barostat
--------------------

DL_POLY_4uses the Melchionna modification of the Nosé-Hoover algorithm
:cite:`melchionna-93a` in which the equations of motion
involve a Nosé-Hoover :index:`thermostat<thermostat;Nosé-Hoover>` 
and a :index:`barostat<barostat;Nosé-Hoover>` in the same spirit.
Additionally, as shown in :cite:`martyna-94a`, a
modification allowing for coupling between the thermostat and barostat
is also introduced.

Cell size variation
~~~~~~~~~~~~~~~~~~~

For isotropic fluctuations the equations of motion are:

.. math::
   :label: npth_eq

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \eta (t) \; (\underline{r}(t) - \underline{R}_{0}(t)) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t)}{m} - \left[ \chi(t) + \eta (t) \right] \underline{v}(t) \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\eta (t)^{2} - 2 \sigma - k_{B}~T_{\rm ext}}{q_{mass}} \nonumber \\
   q_{mass} =& 2~\sigma~\tau_{T}^{2} \label{npth} \\
   \frac{d}{dt}\eta (t) =& 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} - \chi(t) \eta (t) \nonumber \\
   p_{mass} =& (f+3)~k_{B}~T_{\rm ext}~\tau_{P}^{2} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \eta (t) \; \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& [3 \eta (t)]~V(t)~~, \nonumber\end{aligned}

where :math:`\eta` is the barostat friction coefficient,
:math:`\underline{R}_{0}(t)` the system centre of mass at time :math:`t`,
:math:`q_{mass}` the thermostat mass, :math:`\tau_{T}` a specified time
constant for temperature fluctuations, :math:`\sigma` the target
thermostat energy (equation :eq:`sigma_eq`), :math:`p_{mass}`
the barostat mass, :math:`\tau_{P}` a specified time constant for
pressure fluctuations, :math:`{\cal P}` the instantaneous :index:`pressure<units;pressure>`
(equation :eq:`prs_inst_eq`) and :math:`V` the system volume.
is the cell matrix whose columns are the three cell vectors
:math:`\underline{a}, \underline{b}, \underline{c}`.

The conserved quantity is, to within a constant, the Gibbs free energy
of the system:

.. math::

   \begin{aligned}
   {\cal H}_{\rm NPT} = &{\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +{p_{mass}~\eta (t)^{2} \over 2} \\
   & + P_{\rm ext} V(t) +
   (f+1)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~,
   \end{aligned}

where :math:`f` is the system’s degrees of freedom - equation
:eq:`freedom_eq`.

The VV implementation of the Nosé-Hoover algorithm only requires
iterations if bond or PMF constraints are present (:math:`5` until
satisfactory convergence of the constraint forces is achieved). These
are with respect to the pressure (i.e. :math:`\eta (t)`) in the first
part, VV1+RATTLE_VV1. The second part is conventional, VV2+RATTLE_VV2,
as at the end the velocities are scaled by a factor of :math:`\chi`.

#. Thermostat: Note :math:`E_{kin}(t)` changes inside

   .. math::

      \begin{aligned}
      \chi (t + {1 \over 8} \Delta t) \leftarrow& \chi (t) + {\Delta t \over 8} \;
      {{2 E_{kin}(t) + p_{mass}~\eta (t)^{2} - 2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
      \underline{v}(t) \leftarrow& \exp \left( -\chi (t + {1 \over 8} \Delta t) \;
      {\Delta t \over 4} \right) \; \underline{v}(t) \\
      \chi (t + {1 \over 4} \Delta t) \leftarrow& \chi (t + {1 \over 8} \Delta t) + {\Delta t \over 8} \;
      {{2 E_{kin}(t) + p_{mass}~\eta (t)^{2} - 2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber\end{aligned}

#. Barostat: Note :math:`E_{kin}(t)` and :math:`{\cal P}(t)` have
   changed and change inside

   .. math::

      \begin{aligned}
      \eta (t) \leftarrow& \exp \left( -\chi (t + {1 \over 4} \Delta t) \; {\Delta t \over 8} \right) \;
      \eta (t) \nonumber \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \eta (t) + {\Delta t \over 4} \;
      {3 \left[ {\cal P}(t) - P_{\rm ext} \right] V(t) \over p_{mass}} \nonumber \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \exp \left( -\chi (t + {1 \over 4} \Delta t) \;
      {\Delta t \over 8} \right)  \; \eta (t + {1 \over 4} \Delta t) \nonumber \\
      \underline{v}(t) \leftarrow& \exp \left[ -\eta (t + {1 \over 4} \Delta t) \;
      {\Delta t \over 2} \right] \; \underline{v}(t) \\
      \eta (t + {1 \over 4} \Delta t) \leftarrow& \exp \left( -\chi (t + {1 \over 4} \Delta t) \;
      {\Delta t \over 8} \right)  \; \eta (t + {1 \over 4} \Delta t) \nonumber \\
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \eta (t + {1 \over 4} \Delta t) + {\Delta t \over 4} \;
      {3 \left[ {\cal P}(t) - P_{\rm ext} \right] V(t) \over p_{mass}} \nonumber \\
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \exp \left( -\chi (t + {1 \over 4} \Delta t) \;
      {\Delta t \over 8} \right)  \; \eta (t + {1 \over 2} \Delta t) \nonumber\end{aligned}

#. Thermostat: Note :math:`E_{kin}(t)` has changed and changes inside

   .. math::

      \begin{aligned}
      \chi (t + {3 \over 8} \Delta t) \leftarrow& \chi (t + {1 \over 4} \Delta t) + {\Delta t \over 8} \;
      {{2 E_{kin}(t) + p_{mass}~\eta (t + {1 \over 2} \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
      \underline{v}(t) \leftarrow& \exp \left( -\chi (t + {3 \over 8} \Delta t) \;
      {\Delta t \over 4} \right) \; \underline{v}(t) \\
      \chi (t + {1 \over 2} \Delta t) \leftarrow& \chi (t + {3 \over 8} \Delta t) + {\Delta t \over 8} \;
      {{2 E_{kin}(t) + p_{mass}~\eta (t + {1 \over 2} \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber\end{aligned}

#. VV1:

   .. math::

      \begin{aligned}
      \underline{v}(t + {1 \over 2} \Delta t) \leftarrow& \underline{v}(t) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m} \nonumber \\
      \underline{\underline{\mathbf{H}}}(t + \Delta t) \leftarrow& \exp \left[
      \eta (t + {1 \over 2} \Delta t) \; \Delta t \right] \; \underline{\underline{\mathbf{H}}}(t) \nonumber \\
      V(t + \Delta t) \leftarrow& \exp \left[3 \eta (t + {1 \over 2} \Delta t) \;
      \Delta t \right] \; V(t) \\
      \underline{r}(t + \Delta t) \leftarrow& \exp \left[ \eta (t + {1 \over 2} \Delta t) \; \Delta t \right] \;
      (\underline{r}(t) - \underline{R}_{0}(t)) + \Delta t \; \underline{v}(t + {1 \over 2} \Delta t) + \underline{R}_{0}(t) \nonumber\end{aligned}

#. RATTLE_VV1

#. FF:

   .. math:: \underline{f}(t + \Delta t) \leftarrow \underline{f}(t)

#. VV2:

   .. math::

      \begin{aligned}
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + {\Delta t \over 2}) +
      {\Delta t \over 2} \; {\underline{f}(t) \over m}\end{aligned}

#. RATTLE_VV2

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` has changed and
   changes inside

   .. math::

      \begin{aligned}
      \chi (t + {5 \over 8} \Delta t) \leftarrow& \chi (t + {1 \over 2} \Delta t) +
      {\Delta t \over 8} \; {{2 E_{kin}(t + \Delta t) + p_{mass}~\eta (t + {1 \over 2} \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \exp \left(-\chi (t + {5 \over 8} \Delta t) \;
      {\Delta t \over 4} \right) \; \underline{v}(t + \Delta t) \\
      \chi (t + {3 \over 4} \Delta t) \leftarrow& \chi (t + {5 \over 8} \Delta t) +
      {\Delta t \over 8} \; {{2 E_{kin}(t + \Delta t) + p_{mass}~\eta (t + {1 \over 2} \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber\end{aligned}

#. Barostat: Note :math:`E_{kin}(t + \Delta t)` and
   :math:`{\cal P}(t + \Delta t)` have changed and change inside

   .. math::

      \begin{aligned}
      \eta (t + {1 \over 2} \Delta t) \leftarrow& \exp \left( -\chi (t + {3 \over 4} \Delta t) \;
      {\Delta t \over 8} \right) \; \eta (t + {1 \over 2} \Delta t) \nonumber \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \eta (t + {1 \over 2} \Delta t) + {\Delta t \over 4} \;
      {3 \left[ {\cal P}(t + \Delta t) - P_{\rm ext} \right] V(t + \Delta t) \over p_{mass}} \nonumber \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \exp \left( -\chi (t + {3 \over 4} \Delta t) \;
      {\Delta t \over 8} \right) \; \eta (t + {3 \over 4} \Delta t) \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \exp \left[ -\eta (t + {3 \over 4} \Delta t) \;
      {\Delta t \over 2} \right] \; \underline{v}(t + \Delta t) \\
      \eta (t + {3 \over 4} \Delta t) \leftarrow& \exp \left( -\chi (t + {3 \over 4} \Delta t) \;
      {\Delta t \over 8} \right) \; \eta (t + {3 \over 4} \Delta t) \nonumber \\
      \eta (t + \Delta t) \leftarrow& \eta (t + {3 \over 4} \Delta t) + {\Delta t \over 4} \;
      {3 \left[ {\cal P}(t + \Delta t) - P_{\rm ext} \right] V(t + \Delta t) \over p_{mass}} \nonumber \\
      \eta (t + \Delta t) \leftarrow& \exp \left( -\chi (t + {3 \over 4} \Delta t) \;
      {\Delta t \over 8} \right) \; \eta (t + \Delta t) \nonumber\end{aligned}

#. Thermostat: Note :math:`E_{kin}(t + \Delta t)` has changed and
   changes inside

   .. math::

      \begin{aligned}
      \chi (t + {7 \over 8} \Delta t) \leftarrow& \chi (t + {3 \over 4} \Delta t) +
      {\Delta t \over 8} \; {{2 E_{kin}(t + \Delta t) + p_{mass}~\eta (t + \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \exp \left(-\chi (t + {7 \over 8} \Delta t) \;
      {\Delta t \over 4} \right) \; \underline{v}(t + \Delta t) \\
      \chi (t + \Delta t) \leftarrow& \chi (t + {7 \over 8} \Delta t) +
      {\Delta t \over 8} \; {{2 E_{kin}(t + \Delta t) + p_{mass}~\eta (t + \Delta t)^{2} -
      2 \sigma - k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
      \underline{v}(t + \Delta t) \leftarrow& \underline{v}(t + \Delta t) - \underline{V}_{0}(t + \Delta t)~~, \nonumber\end{aligned}

where :math:`\underline{V}_{0}(t + \Delta t)` is the c.o.m. velocity at
timestep :math:`t + \Delta t` and is the cell matrix whose columns are
the three cell vectors :math:`\underline{a}, \underline{b}, \underline{c}`.

The VV flavour of the Nosé-Hoover barostat (and thermostat) is
implemented in the DL_POLY_4routine ``npt_h0_vv``. The routine
``npt_h1_vv`` implements the same but also incorporate RB dynamics.

Cell size and shape variations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The isotropic algorithmscmay be extended to allowing the cell shape to
vary by defining :math:`\eta` as a tensor, :math:`\underline{\underline{\mathbf{\eta}}}`. The
equations of motion are written in the same fashion as is in the
isotropic algorithm with slight modifications (as now the equations with
:math:`\eta` are extended to matrix forms)

.. math::
   :label: nsth_eq

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \underline{\underline{\mathbf{\eta}}}(t) \cdot (\underline{r}(t) - \underline{R}_{0}(t)) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t)}{m} - \left[ \chi(t)~\underline{\underline{\mathbf{1}}} +
   \underline{\underline{\mathbf{\eta}}}(t) \right] \cdot \underline{v}(t) \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - 3^{2}~k_{B}~T_{\rm ext}}{q_{mass}} \nonumber \\
   q_{mass} =& 2~\sigma~\tau_{T}^{2} \label{nsth} \\
   \frac{d}{dt}\underline{\underline{\mathbf{\eta}}}(t) =& \frac{\underline{\underline{\mathbf{\sigma}}}(t) -
   P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}} - \chi(t) \underline{\underline{\mathbf{\eta}}}(t) \nonumber \\
   p_{mass} =& \frac{(f+3)}{3}~k_{B}~T_{\rm ext}~\tau_{P}^{2} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \underline{\underline{\mathbf{\eta}}}(t) \cdot \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& \texttt{ Tr} [\underline{\underline{\mathbf{\eta}}}(t)]~V(t)~~, \nonumber\end{aligned}

where :math:`\underline{\underline{\mathbf{\sigma}}}` is the stress tensor
(equation :eq:`str_inst_eq`) and :math:`\underline{\underline{\mathbf{1}}}` is the
identity matrix. The VV algorithmic equations are, therefore, written in
the same fashion as above with slight modifications in (i) the equations
for the thermostat and barostat frictions, and (ii) the equations for
the system volume and cell parameters. The modifications in (i) for the
VV couched algorithm are of the following sort

.. math::

   \begin{aligned}
   \chi (t + {1 \over 8} \Delta t) \leftarrow& \chi (t) + {\Delta t \over 8} \;
   {{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - 3^{2}~k_{B}~T_{\rm ext}} \over q_{mass}} \nonumber \\
   \underline{v}(t) \leftarrow& \exp \left[ - \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4}\Delta t) \;
   {\Delta t \over 2} \right] \cdot \underline{v}(t) \\
   \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4} \Delta t) \leftarrow& \underline{\underline{\mathbf{\eta}}}(t) +
   {\Delta t \over 4} \; \frac{\underline{\underline{\mathbf{\sigma}}}(t) - P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}}~~, \nonumber\end{aligned}

The modifications in (ii) couched algorithms

.. math::

   \begin{aligned}
   \underline{\underline{\mathbf{H}}}(t + \Delta t) \leftarrow& \exp \left( \underline{\underline{\mathbf{\eta}}}(t + {1 \over 2} \Delta t) \;
   \Delta t \right) \cdot \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   V(t + \Delta t) \leftarrow& \exp \left(\texttt{ Tr}
   \left[ \underline{\underline{\mathbf{\eta}}}(t + {1 \over 2} \Delta t) \right] \; \Delta t \right) \; V(t)~~.\end{aligned}

It is worth noting DL_POLY_4uses Taylor expansion truncated to the
quadratic term to approximate exponentials of tensorial terms.

The conserved quantity is, to within a constant, the Gibbs free energy
of the system:

.. math::

   \begin{aligned}
   {\cal H}_{\rm N\underline{\underline{\mathbf{\sigma}}}T} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2}  \\ 
   &~~~+P_{\rm ext} V(t) +
   (f+3^{2})~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~,
   \end{aligned}

where :math:`f` is the system’s degrees of freedom - equation
:eq:`freedom_eq`.

This ensemble is optionally extending to constant normal pressure and
constant surface area, NP\ :math:`_{n}`\ AT
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion and slight amending the thermostat equation
of motion and the conserved quantity to:

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\beta}(t) =& \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} -
   \chi(t) \eta_{zz}(t) & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha,\beta) \ne z
   \end{array} \right. \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - k_{B}~T_{\rm ext}}{q_{mass}} \\
   {\cal H}_{\rm NP_{n}AT} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\
   &~~+ P_{\rm ext} V(t) +
   (f+1)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~. \nonumber\end{aligned}

Similarly, this ensemble is optionally extending to constant normal
pressure and constant surface tension, NP\ :math:`_{n}\gamma`\ T
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion and slight amending the thermostat equation
of motion and the conserved quantity to:

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\beta}(t) =& \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{\alpha\alpha}(t) - \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} -
   \chi(t) \eta_{\alpha\alpha}(t) & (\alpha = \beta) = x,y \\
   & \\
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} -
   \chi(t) \eta_{zz}(t) & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha \ne \beta) = x,y,z
   \end{array} \right. \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - 3~k_{B}~T_{\rm ext}}{q_{mass}} \\
   {\cal H}_{\rm NP_{n}\gamma T} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\
   &~~~+ P_{\rm ext} V(t) +
   (f+3)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~. \nonumber\end{aligned}

where :math:`\gamma_{\rm ext}` is the user defined external surface
tension and :math:`h_{z}(t) = V(t) / A_{xy}(t)` is the instantaneous
hight of the MD box (or MD box volume over area). One defines the
instantaneous surface tension as given in
equation :eq:`gamma_eq`. The case :math:`\gamma_{\rm ext}=0`
generates the NPT anisotropic ensemble for the orthorhombic cell
(\ ``imcon``\ =2 in CONFIG, see
:ref:`Appendix B<boundary-conditions>`). This can
be considered as an "orthorhombic" constraint on the
N\ :math:`\sigma`T ensemble. The constraint can be strengthened
further, to a "semi-orthorhombic" one, by imposing that the MD cell
change isotropically in the :math:`(x,y)` plane which leads to the
following changes in the equations above

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\alpha}(t) =& \frac{\left[\sigma_{xx}(t)+\sigma_{yy}(t)\right]/2 -
   \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} -
   \chi (t) \eta_{\alpha\alpha}(t)~~:~~(\alpha = \beta) = x,y \nonumber \\
    & \\
   {\cal H}_{\rm NP_{n}{\gamma}{=0}T} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\ 
   &~~~+ P_{\rm ext} V(t) +
   (f+2)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~.\nonumber\end{aligned}

The VV flavour of the non-isotropic Nosé-Hoover barostat (and
thermostat) is implemented in the DL_POLY_4routine ``nst_h0_vv``. The
routine ``nst_h1_vv`` implements the same but also incorporate RB
dynamics.

Martyna-Tuckerman-Klein Barostat
--------------------------------

DL_POLY_4includes the Martyna-Tuckerman-Klein (MTK) interpretation of
the VV flavoured Nosé-Hoover algorithms :cite:`martyna-96a`
for isotropic and anisotropic cell fluctuations in which the equations
of motion are only slightly augmented with respect to those for the
coupled Nosé-Hoover :index:`thermostat<thermostat;Nosé-Hoover>` and 
:index:`barostat<barostat;Nosé-Hoover>`. Compare the isotropic cell
changes case, equations :eq:`npth_eq`, to

.. math::

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \eta (t) \; \underline{r}(t) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t)}{m} - \left[ \chi(t) +
   \left(1+\frac{3}{f}\right) \eta (t) \right] \underline{v}(t) \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\eta (t)^{2} -
   2 \sigma - k_{B}~T_{\rm ext}}{q_{mass}} \nonumber \\
   q_{mass} =& 2~\sigma~\tau_{T}^{2} \\
   \frac{d}{dt}\eta (t) =& 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} +
   3 \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} - \chi(t) \eta (t) \nonumber \\
   p_{mass} =& (f+3)~k_{B}~T_{\rm ext}~\tau_{P}^{2} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \eta (t) \; \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& [3 \eta (t)]~V(t)~~, \nonumber\end{aligned}

and the anisotropic cell change case, equation :eq:`nsth_eq`,
to

.. math::

   \begin{aligned}
   \frac{d}{dt} \underline{r}(t) =& \underline{v}(t) + \underline{\underline{\mathbf{\eta}}}(t) \cdot  \underline{r}(t) \nonumber \\
   \frac{d}{dt} \underline{v}(t) =& \frac{\underline{f}(t)}{m} - \left[ \chi(t)~\underline{\underline{\mathbf{1}}} +
   \underline{\underline{\mathbf{\eta}}}(t) + \frac{\texttt{ Tr}\left[\underline{\underline{\mathbf{\eta}}}(t)\right]}{f}~\underline{\underline{\mathbf{1}}} \right] \cdot \underline{v}(t) \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - 3^{2}~k_{B}~T_{\rm ext}}{q_{mass}} \nonumber \\
   q_{mass} =& 2~\sigma~\tau_{T}^{2} \\
   \frac{d}{dt}\underline{\underline{\mathbf{\eta}}}(t) =& \frac{\underline{\underline{\mathbf{\sigma}}}(t) -
   P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}} + \frac{2 E_{kin}(t)}{f} \frac{\underline{\underline{\mathbf{1}}}}{p_{mass}} -
   \chi(t) \underline{\underline{\mathbf{\eta}}}(t) \nonumber \\
   p_{mass} =& \frac{(f+3)}{3}~k_{B}~T_{\rm ext}~\tau_{P}^{2} \nonumber \\
   \frac{d}{dt}\underline{\underline{\mathbf{H}}}(t) =& \underline{\underline{\mathbf{\eta}}}(t) \cdot \underline{\underline{\mathbf{H}}}(t) \nonumber \\
   \frac{d}{dt} V(t) =& \texttt{ Tr} [\underline{\underline{\mathbf{\eta}}}(t)]~V(t)~~. \nonumber\end{aligned}

The changes include one extra dependence to the velocity and barostat
equations and removal of the centre of mass variable
:math:`\underline{R}_{0}(t)` dependence in the position equation.

The modifications in for the VV couched algorithms are of the following
sort

.. math::

   \begin{aligned}
   \eta (t + {1 \over 4} \Delta t) \leftarrow& \eta (t) + {\Delta t \over 4} \;
   \left[ 3 V(t) \frac{{\cal P}(t) - P_{\rm ext}}{p_{mass}} +
   3 \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} \right] \nonumber \\
   \underline{v}(t) \leftarrow& \exp \left[ -\left( 1 + \frac{3}{f} \right)
   \eta (t + {1 \over 4}\Delta t) \; {\Delta t \over 2} \right] \; \underline{v}(t) \\
   \underline{r}(t + \Delta t) \leftarrow& \exp \left[ \eta (t + {1 \over 2} \Delta t) \; \Delta t \right] \;
   \underline{r}(t) + \Delta t \; \underline{v}(t + {1 \over 2} \Delta t) \nonumber\end{aligned}

for the isotropic cell fluctuations case and

.. math::

   \begin{aligned}
   \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4} \Delta t) \leftarrow& \underline{\underline{\mathbf{\eta}}}(t) +
   {\Delta t \over 4} \; \left[ \frac{\underline{\underline{\mathbf{\sigma}}}(t) - P_{\rm ext}~V(t)~\underline{\underline{\mathbf{1}}}}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f} \frac{\underline{\underline{\mathbf{1}}}}{p_{mass}} \right] \nonumber \\
   \underline{v}(t) \leftarrow& \exp \left[ -\left( \underline{\underline{\mathbf{\eta}}}(t + {1 \over 4}\Delta t) +
   \frac{1}{f} \texttt{ Tr}\left[\underline{\underline{\mathbf{\eta}}}(t + {1 \over 4}\Delta t)\right] \right) \;
   {\Delta t \over 2} \right] \cdot \underline{v}(t) \\
   \underline{r}(t + \Delta t) \leftarrow& \exp \left[ \underline{\underline{\mathbf{\eta}}} (t + {1 \over 2} \Delta t) \; \Delta t \right] \cdot
   \underline{r}(t) + \Delta t \; \underline{v}(t + {1 \over 2} \Delta t) \nonumber\end{aligned}

for the anisotropic cell fluctuations case.

This ensemble is optionally extending to constant normal pressure and
constant surface area, NP\ :math:`_{n}`\ AT
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion and slight amending the thermostat equation
of motion and the conserved quantity to:

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\beta}(t) =& \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} + \frac{2 E_{kin}(t)}{f~p_{mass}} -
   \chi(t) \eta_{zz}(t) & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha,\beta) \ne z
   \end{array} \right. \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - k_{B}~T_{\rm ext}}{q_{mass}} \\
   {\cal H}_{\rm NP_{n}AT} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\
   &~~~+ P_{\rm ext} V(t) +
   (f+1)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~. \nonumber\end{aligned}

Similarly, this ensemble is optionally extending to constant normal
pressure and constant surface tension, NP\ :math:`_{n}\gamma`\ T
:cite:`ikeguchi-04a`, by semi-isotropic constraining of the
barostat equation of motion and slight amending the thermostat equation
of motion and the conserved quantity to:

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\beta}(t) =& \left\{ \begin{array} {l@{\quad:\quad}l}
   \frac{\sigma_{\alpha\alpha}(t) - \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f} \frac{1}{p_{mass}} - \chi(t) \eta_{\alpha\alpha}(t) & (\alpha = \beta) = x,y \\
   & \\
   \frac{\sigma_{zz}(t) - P_{\rm ext}~V(t)}{p_{mass}} +
   \frac{2 E_{kin}(t)}{f~p_{mass}} - \chi(t) \eta_{zz}(t) & (\alpha = \beta) = z \\
   0~~~;~~~\eta_{\alpha\beta}(0) = 0 & (\alpha \ne \beta) = x,y,z
   \end{array} \right. \nonumber \\
   \frac{d}{dt} \chi(t) =& \frac{2 E_{kin}(t) + p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}}(t) \cdot
   \underline{\underline{\mathbf{\eta}}}(t)^{T}] - 2 \sigma - 3~k_{B}~T_{\rm ext}}{q_{mass}} \\
   {\cal H}_{\rm NP_{n}\gamma T} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\
   &~~~+ P_{\rm ext} V(t) +
   (f+3)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~, \nonumber\end{aligned}

where :math:`\gamma_{\rm ext}` is the user defined external surface
tension and :math:`h_{z}(t) = V(t) / A_{xy}(t)` is the instantaneous
hight of the MD box (or MD box volume over area). One defines the
instantaneous surface tension as given in
equation :eq:`gamma_eq`. The case :math:`\gamma_{\rm ext}=0`
generates the NPT anisotropic ensemble for the orthorhombic cell
(\ ``imcon``\ =2 in CONFIG, see
:ref:`Appendix B<boundary-conditions>`). This can
be considered as an "orthorhombic" constraint on the
N\ :math:`\sigma`T ensemble. The constraint can be strengthened
further, to a "semi-orthorhombic" one, by imposing that the MD cell
change isotropically in the :math:`(x,y)` plane which leads to the
following changes in the equations above

.. math::

   \begin{aligned}
   \frac{d}{dt} \eta_{\alpha\alpha}(t) =& \frac{\left[\sigma_{xx}(t)+\sigma_{yy}(t)\right]/2 -
   \left[ P_{\rm ext} - \gamma_{\rm ext} / h_{z}(t) \right]~V(t)}{p_{mass}} + \frac{2~E_{kin}(t)}{f~p_{mass}} -
   \chi (t) \eta_{\alpha\alpha}(t)~:~\nonumber \\
    & \phantom{xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx}~:~(\alpha = \beta) = x,y\\
   {\cal H}_{\rm NP_{n}{\gamma}{=0}T} =& {\cal H}_{\rm NVE} + {q_{mass}~\chi (t)^{2} \over 2} +
   {p_{mass}~\texttt{ Tr}[\underline{\underline{\mathbf{\eta}}} \cdot \underline{\underline{\mathbf{\eta}}}^{T}] \over 2} \\
   &~~~+ P_{\rm ext} V(t) +
   (f+2)~k_{B}~T_{\rm ext}~\int_o^t \chi (s) ds~~.\nonumber\end{aligned}

Although the Martyna-Tuckerman-Klein equations of motion have same
conserved quantities as the Nosé-Hoover’s ones they are proven to
generate ensembles that conserve the phase space volume and thus have
well defined conserved quantities even in presence of forces external to
the system :cite:`martyna-94a`, which is not the case for
Nosé-Hoover NPT and N\ :math:`\underline{\underline{\mathbf{\sigma}}}`\ T ensembles.

The NPT and N\ :math:`\underline{\underline{\mathbf{\sigma}}}`\ T versions of the MTK ensemble are
implemented in the DL_POLY_4routines ``npt_m0_vv`` and ``nst_m0_vv``.
The corresponding routines incorporating RB dynamics are ``npt_m1_vv``,
and ``nst_m1_vv``.
