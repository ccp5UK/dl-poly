.. _DPD_barostats:

Barostats 
=========

In addition to a DPD thermostat, a DPD based barostat may be included in simulations to obtain a desired 
average pressure (:math:`P_{0}`) by adjusting the size (and shape) of the simulation cell. These barostats 
are related to the typical molecular dynamics versions adapted for the DPD formalism. 

The isotropic pressure in a system is calculated using the virial theorem, 

.. math::

    P(t) = \frac{1}{3V(t)} \left[ \sum_{i}m_{i}v_{i}^{2}\left(t\right) + \sum_{i}\underline{F}_{i}\left(t\right) \cdot \underline{r}_{i}\left(t\right) \right],

while for anisotropic orthorhombic systems the pressure in dimension :math:`\alpha`, equal to the 
instantaneous stress tensor component :math:`\sigma_{\alpha\alpha}\left(t\right)^{2}`, is identified as, 

.. math::

    P_{\alpha}\left(t\right) = \frac{1}{V\left(t\right)} \left[ \sum_{i}m_{i}v_{i,\alpha}^{2}\left(t\right) + \sum_{i}F_{i,\alpha}\left(t\right)r_{i,\alpha}\left(t\right) \right] = \sigma_{\alpha\alpha}.

In both equations, the instantaneous values required for barostats include only the interaction forces 
(e.g. soft pairwise interactions, bonds, electrostatics), they do not include virial contributions from 
thermostatting, which are included in reported values of system pressure. 

.. note::
    The stress tensor :math:`\sigma_{\alpha\beta}` is equal at least in magnitude to the pressure tensor 
    :math:`P_{\alpha\beta}`\ : the two tensors may be equal to each other depending on sign convention. 
    For analysis purposes, the two tensors may be used interchangeably. 


Langevin Barostat 
-----------------

The governing equation for the Langevin barostat on an orthorhombic simulation cell :cite:`quigley2004,jakobsen2005a,jakobsen2005aerratum` is 
the force exerted by the piston (expressed as the time-derivative of its momentum 
:math:`p_{g,\alpha} = W_{g}u_{g,\alpha}`):

.. math::
    :label: DPD_langevin_barostat_eq 

    \dot{p}_{g,\alpha} &= V \left( P_{\alpha} - P_0 + \frac{\gamma_0}{L_z} \right) + \frac{1}{N_f} \sum_i m_i v_i^2 - \gamma_p p_{g,\alpha} + \sigma_p \zeta_{p, \alpha} \Delta t^{-\frac{1}{2}} ~~~ (\alpha = x, y) \\
    \dot{p}_{g,z} &= V \left( P_{z} - P_0 \right) + \frac{1}{N_f} \sum_i m_i v_i^2 - \gamma_p p_{g,z} + \sigma_p \zeta_{p,g} \Delta t^{-\frac{1}{2}}

where :math:`N_{f}` is the number of degrees of freedom, i.e. for a three-dimensional box containing :math:`N` 
moving particles :math:`N_{f} = 3(N - 1)`. :math:`\gamma_{p}` and :math:`\sigma_{p} \equiv \sqrt{\tfrac{2}{3} \gamma_p W_g k_B T}` 
are respectively the drag and random coefficients for the piston and :math:`\zeta_{p,\alpha}` is a Gaussian 
random number for dimension :math:`\alpha`. All three components of the barostat momentum :math:`\underline{p}_{g}`
and the Gaussian random numbers can be set equal to the same values if operating isotropically, while the 
:math:`x`- and :math:`y`-components for the momentum and corresponding random numbers can be set to the 
same values when operating semi-isotropically, or zero if using a constant surface area ensemble (as only
the :math:`z`-component of the piston momentum needs to evolve in this case). When both :math:`\gamma_{p}` 
and :math:`\sigma_{p}` are set to zero, the Langevin barostat reduces to the **extended system** method. 

The subsequent simulation cell size :math:`L_{\alpha}` can be determined by 

.. math::
    :label: DPD_langevin_cell_size_eq 

    \dot{L}_{\alpha} = \frac{p_{g,\alpha} L_{\alpha}}{W_g} = u_{g,\alpha} L_{\alpha}

with the barostat mass :math:`W_{g}` chosen to be equal to :math:`Nk_{\text{B}}T\tau_{p}^{2}` where 
:math:`\tau_{p}` is the characteristic barostat time and should be set equal to between :math:`\frac{2}{\gamma_{p}}` 
and :math:`\frac{10}{\gamma_{p}}`. 

The velocities and positions of the particles are calculated by integration of slightly modified differential 
equations, 

.. math::
    :label: DPD_langevin_dv_eq 

    \frac{dv_{i,\alpha}}{dt} = \frac{F_{i,alpha}}{m_{i}} - u_{g,\alpha}v_{i,\alpha} - \frac{1}{N_{f}}v_{i,\alpha} \sum_{\alpha}u_{g,\alpha},

.. math::
    :label: DPD_langevin_dr_eq 
    
    \frac{dr_{i,\alpha}}{dt} = v_{i,\alpha} + u_{g,\alpha}r_{i,\alpha}

where the force on particle :math:`i`, :math:`\underline{F}_{i}`, includes any thermostatting forces: the 
time integral of these forces can be determined for all thermostat types. 

The implementation of this barostat is dependent on being able to integrate the equation of motion for the 
barostat alongside those for the particles. It is possible to do this by using the standard Velocity Verlet 
scheme to integrate the barostat velocity, :cite:`jakobsen2005a,jakobsen2005aerratum` although the implicit dependence of barostat 
force on the barostat velocity requires iteration of these properties and the particle velocities during 
the second stage to ensure self-consistency. To avoid iteration and enable to use the barostat in 
circumstances when convergence cannot be achieved, it is possible to devise an alternate explicit method 
of integrating the barostat velocity using Trotter factorisation :cite:`martyna1996,quigley2004`. 

.. note::
    The modification of this approach to incorporate pairwise thermostats (e.g. DPD) as shown here was 
    originally carried out by Yaser Afshar, Freiderike Schmid and Mohammad Said Saidi in 2013 as provided 
    in a hitherto unpublished manuscript ('An explicit derivation of Langevin piston approach for constant-
    pressure dissipative particle dynamics').

The following Liouville propagation operator for the equations of motion including the Langevin barostat can 
be expressed as, 

.. math::
    
    \mathcal{L}= \mathcal{L}_{r} + \mathcal{L}_{P} + \mathcal{L}_{\epsilon_{\alpha}} + \mathcal{L}_{NHP_{g}} 

where the individual operators for a given dimension of the piston (:math:`\alpha`) are given as 

.. math::

    \mathcal{L}_{r} &= \sum_{i}\left[\underline{v}_{i} + \frac{p_{g,\alpha}}{W_{g}}\underline{r}_{i}\right] \cdot \nabla_{\underline{r}_{i}} \\
    \mathcal{L}_{P} &= \sum_{i} \left[ \frac{\underline{F}_{i}}{m_{i}} \right] \cdot \nabla_{\underline{v}_{i}} \\
    \mathcal{L}_{\epsilon_{\alpha}} &= \frac{p_{g,\alpha}}{W_{g}}\frac{\partial}{\partial\epsilon_{\alpha}} \\ 
    \mathcal{L}_{NHP_{g}} &= -\frac{1}{W_{g}} \left( p_{g,\alpha} + \frac{\sum_{\alpha}p_{g,\alpha}}{N_{f}} \right) \sum_{i} \underline{v}_{i}\cdot\nabla_{\underline{v}_{i}} + \left[G_{\epsilon} - \gamma_{p}p_{g,\alpha}\right] \frac{\partial}{\partial p_{g,\alpha}}

with :math:`\epsilon_{\alpha} = \ln L_{\alpha}` and 

.. math::

    G_{\epsilon} = \left[ \left( 1 + \frac{1}{N_{f}} \right) \sum_{i}m_{i}v_{i}^{2} + \sum_{i}\underline{r}_{i}\cdot\underline{F}_{i} - VP_{0} + \sigma_{p}\zeta_{p,\alpha}\Delta t^{-\frac{1}{2}} \right].

These equation of motion can be integrated using the following Trotter expansion, 

.. math:: 

    \mathcal{e}^{\left(\mathcal{L}\Delta t\right)} = \mathcal{e}^{\left(\mathcal{L}_{NHP_{g}}\frac{\Delta t}{2}\right)} \mathcal{e}^{\left(\mathcal{L}_{P}\frac{\Delta t}{2}\right)} \mathcal{e}^{\left(\mathcal{L}_{\epsilon_{\alpha}}\Delta t\right)} \mathcal{e}^{\left(\mathcal{L}_{r}\Delta t\right)} \mathcal{e}^{\left(\mathcal{L}_{P}\frac{\Delta t}{2}\right)} \mathcal{e}^{\left(\mathcal{L}_{NHP_{g}}\frac{\Delta t}{2}\right)}

The part with :math:`\mathcal{L}_{NHP_{g}}` can be simplified into, :cite:`martyna1996`

.. math::

    \begin{aligned}
    \mathcal{e}^{\left(\mathcal{L}_{NHP_{g}}\frac{\Delta t}{2}\right)} &= \mathcal{e}^{\left( -\frac{\Delta t}{8}\gamma_{p}p_{g,\alpha}\frac{\partial}{\partial p_{g,\alpha}} \right)} \mathcal{e}^{\left( \frac{\Delta t}{4}G_{\epsilon}\frac{\partial}{\partial p_{g,\alpha}} \right)} \mathcal{e}^{\left( -\frac{\Delta t}{8}\gamma_{p}p_{g,\alpha}\frac{\partial}{\partial p_{g,\alpha}} \right)} \\ 
    & \times \mathcal{e}^{\left( -\frac{\Delta t}{2}\left(1 + 3N_{f}\right) \frac{p_{g,\alpha}}{W_{g}}\underline{v}_{i}\cdot\nabla_{\underline{v}_{i}}\right)} \\
    & \times \mathcal{e}^{\left( -\frac{\Delta t}{8}\gamma_{p}p_{g,\alpha}\frac{\partial}{\partial p_{g,\alpha}} \right)} \mathcal{e}^{\left( \frac{\Delta t}{4}G_{\epsilon}\frac{\partial}{\partial p_{g,\alpha}} \right)} \mathcal{e}^{\left( -\frac{\Delta t}{8}\gamma_{p}p_{g,\alpha}\frac{\partial}{\partial p_{g,\alpha}} \right)}
    \end{aligned}

Each of the operators shown above can be expressed as a simple calculation using the most recent values of 
piston and particle velocities. 

The first stage of integration starts by advancing the piston velocities by quarter of a timestep: in line with
the Trotter expansion for the barostat itself (:math:`\mathcal{L}_{NHP_{g}}`), this can be carried out in three
stages, 

.. math::

    p_{g,\alpha}\left(t^{+}\right) = \mathcal{e}^{-\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha}\left(r\right) 


.. math::

    p_{g,\alpha}\left(t^{++}\right) = p_{g,\alpha}\left(t^{+}\right) + \frac{\Delta t}{4} \left[ V\left(t\right)\left(P_{\alpha}\left(t\right) - P_{0} + \frac{\gamma_{0}}{L_{z}}\right) + \frac{1}{N_{f}} \sum_{i}m_{i}v{i}^{2}\left(t\right) + \sigma_{p}\zeta_{p,\alpha}\left(t\right)\Delta t^{\frac{1}{2}} \right] 


.. math::

    p_{g,\alpha}\left(t + \frac{1}{4}\Delta t \right) = \mathcal{e}^{\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha} \left(t^{++}\right)


The particle velocities are then shifted to the half timestep by integrating forces in the usual manner, 

.. math::

    \underline{v}_{i}\left(t + \frac{t}{2}\Delta t\right) = \underline{v}_{i}\left(t^{+}\right) + \frac{\Delta t}{2 m_{i}}\underline{F}_{i}\left(t\right)

The piston velocity is then advanced by another quater timestep, using the most recently updated particle 
velocities, 

.. math::

    p_{g,\alpha}\left(t^{+} + \frac{1}{4}\Delta t\right) = \mathcal{e}^{-\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha}\left(r\right) 


.. math::

    p_{g,\alpha}\left(t^{++} + \frac{1}{4}\Delta t\right) = p_{g,\alpha}\left(t^{+} + \frac{1}{4}\Delta t\right) + \frac{\Delta t}{4} \left[ V\left(t\right)\left(P_{\alpha}\left(t\right) - P_{0} + \frac{\gamma_{0}}{L_{z}}\right) + \frac{1}{N_{f}} \sum_{i}m_{i}v{i}^{2}\left(t\right) + \sigma_{p}\zeta_{p,\alpha}\left(t\right)\Delta t^{\frac{1}{2}} \right] 


.. math::

    p_{g,\alpha}\left(t + \frac{1}{2}\Delta t \right) = \mathcal{e}^{\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha} \left(t^{++} + \frac{1}{4}\Delta t\right)

while the particle positions are simultaneously rescaled due to the change in box volume and advanced to the 
end of the timestep using the misstep velocities, 

.. math::

    r_{i,\alpha}\left(t + \Delta t\right) = r_{i,\alpha}\left(t\right) \exp\left[\frac{p_{g,\alpha}\left(t + \frac{1}{2}\Delta t\right)}{W_{g}}\Delta t\right] + v_{i,\alpha}\left(t + \frac{\Delta t}{2}\right) \Delta t. 

The same exponential factor is used to scale the dimensions of the box, 

.. math:: 

    L_{\alpha}(t + \Delta t) = L_{\alpha}(t) \exp\left[\frac{p_{g,\alpha}\left(t + \frac{1}{2}\Delta t\right)}{W_{g}}\Delta t\right]

and any required boundary conditions can then be applied. The positions of particles involved in fixed-length 
constraints (e.g. applying the first RATTLE stage) can be applied at this point, although the resulting 
additional stress tensor and virial requires an update to the instantaneous pressure. Repeated applications 
of the first integration stage are therefore needed if constraints are included.

After particle forces and the instantaneous pressure for the end of the timestep have been calculated, a new 
random number for the barostat :math:`\zeta_{p,\alpha}\left(t\right)` is generated. The second stage of 
integration starts with integration of the particle forces to give updated particle velocities in the usual 
manner for the second Velocity Verlet stage, 

.. math::

    \underline{v}_{i}\left(t^{+} + \frac{1}{2}\Delta t\right) = \underline{v}_{i}\left(t^{+}\right) + \frac{\Delta t}{2m_{i}}\underline{F}_{i}\left(t + \Delta t\right),

before any require pairwise thermostatting (e.g. second-order Shardlow splitting) is applied at this point, 
along with the second RATTLE stage for any fixed-length constraints. The piston velocity is then advanced 
by another quarter timestep, 

.. math::

    p_{g,\alpha}\left(t^{+} + \frac{1}{2}\Delta t\right) = \mathcal{e}^{-\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha}\left(r\right) 


.. math::

    p_{g,\alpha}\left(t^{++} + \frac{1}{2}\Delta t\right) = p_{g,\alpha}\left(t^{+} + \frac{1}{2}\Delta t\right) + \frac{\Delta t}{4} \left[ V\left(t\right)\left(P_{\alpha}\left(t\right) - P_{0} + \frac{\gamma_{0}}{L_{z}}\right) + \frac{1}{N_{f}} \sum_{i}m_{i}v{i}^{2}\left(t\right) + \sigma_{p}\zeta_{p,\alpha}\left(t\right)\Delta t^{\frac{1}{2}} \right] 


.. math::

    p_{g,\alpha}\left(t + \frac{3}{4}\Delta t \right) = \mathcal{e}^{\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha} \left(t^{++} + \frac{1}{2}\Delta t\right)

before the particle velocities are rescaled again, 

.. math::

    \underline{v}_{i}\left(t^{+} + \Delta t\right) = \underline{v}_{i}\left(t^{+} + \frac{1}{2}\Delta t\right) \exp\left[ - \frac{\Delta t}{2} \left(p_{g,\alpha}\left(t + \frac{3}{4}\Delta t\right) + \frac{\sum_{\alpha}p_{g,\alpha}\left(t + \frac{3}{4}\Delta t\right)}{N_{f}}\right) \right]. 

and the piston velocity is advanced by another quarter timestep, 

.. math::

    p_{g,\alpha}\left(t^{+} + \frac{3}{4}\Delta t\right) = \mathcal{e}^{-\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha}\left(r\right) 


.. math::

    p_{g,\alpha}\left(t^{++} + \frac{3}{4}\Delta t\right) = p_{g,\alpha}\left(t^{+} + \frac{3}{4}\Delta t\right) + \frac{\Delta t}{4} \left[ V\left(t\right)\left(P_{\alpha}\left(t\right) - P_{0} + \frac{\gamma_{0}}{L_{z}}\right) + \frac{1}{N_{f}} \sum_{i}m_{i}v{i}^{2}\left(t\right) + \sigma_{p}\zeta_{p,\alpha}\left(t\right)\Delta t^{\frac{1}{2}} \right] 


.. math::

    p_{g,\alpha}\left(t + \Delta t \right) = \mathcal{e}^{\gamma_{p}\frac{\Delta t}{8}} p_{g,\alpha} \left(t^{++} + \frac{3}{4}\Delta t\right)

.. note::
    No adjustment of particle forces or positions is required during the second integration stage, nor 
    is repeated application of this integration stage. 

To use the DPD based Langevin barostat the ``ensemble`` directive in the ``CONTROL`` file must be set to 
either ``npt`` or ``nst`` and the ``ensemble_method`` must be set to ``dpd``. The ``ensemble_dpd_order`` is 
then used to specify the coupled thermostat as usual. The parameters :math:`\sigma_{p}` and :math:`\gamma_{p}` 
are controlled by the ``ensemble_barostat_coupling`` and ``ensemble_barostat_friction`` directives respectively. 