.. _heat-flux:

Heat Flux
=========

Introduction
~~~~~~~~~~~~

It is possible to use DL_POLY_5 to calculate the heat flux of a material
with two-body interactions in an MD simulation. The heat flux can
subsequently be used as a means to calculate the thermal conductivity of
a material via the Green-Kubo relation (see Section :ref:`correlation-functions`). 

It is also possible to calculate the momentum density, which may be required to apply corrections to the 
thermal-conductivity in multi-component systems when comparing to experimental data.

To enable the calculation of heat flux add **heat_flux On** into the CONTROL file. For the momentum
density supply a list of atom types to calculate the values for e.g. **momentum_density [Li F]** to
compute for Li and F atoms.

Heat flux is currently supported for: direct and tabulated VDW interactions (see Table :numref:`(%s)<vdw-table>`), SPME interactions, and direct and tabulated metal potentials (see Table :numref:`(%s)<metal-table>`). Tersoff, three, and four body potentials are currently unsupported.

Theory
~~~~~~

The heat flux for two-body interactions is defined as:

.. math:: \mathbf{J} = \frac{1}{V} \sum\limits^{N}_{i} \left[ e_{i} \mathbf{v}_{i} + \frac{1}{2}\sum\limits_{j\neq i} \mathbf{f}_{ij}\cdot\mathbf{v}_{i}\mathbf{r}_{ij} \right]
   :label: heat-flux-definition_eq

where :math:`\mathbf{\mathbf{\textbf{J}}}` is the heat flux,
:math:`V` is the volume of the cell, :math:`N` is the number of
particles, :math:`e_{i}` is the energy, :math:`\mathbf{v}_{i}` is the
velocity, :math:`\mathbf{f}_{ij}` is the force of j on i. All
subscripts :math:`i` refer to a particle.

The thermal conductivity can then be calculated as an auto-correlation of the heat
flux

.. math:: \kappa = \frac{V}{k_{B} T^{2}} \int\limits_{0}^{\infty} \langle \mathbf{J}(0)  \mathbf{J}(t) \rangle \, \mathrm{d}t

where :math:`\kappa` is the thermal conductivity, :math:`V` is the
volume of the cell, :math:`k_{B}` is the Boltzmann constant.

In mixed systems where components significantly differ in mass it may be 
necessary to correct for mass transfer effects when calculating 
thermal-conductivity :cite:`armstrong2014thermal`. In order to 
facilitate these corrections the momentum density may also be 
calculated. This is defined for components (atomic types) C as,

.. math:: \mathbf{q}_{C} = \frac{1}{V}\sum_{i\in C}m_{i}\mathbf{v}_{i}.

Where :math:`C` is the set of atoms.

Implementation
~~~~~~~~~~~~~~

For the purposes of calculating per-particle SPME interactions, the
long-range electrostatics forces and energies are calculated differently
for the heat flux case. From the SPME equations, we can calculate a per
particle contribution via:

.. math::

   \begin{gathered}\Omega ^{\mathit{ABC}}=\sum _j\omega _j^{\mathit{ABC}}\\\omega _j^{\mathit{ABC}}=\frac 1{2\pi V}\sum
   _{n_1n_2n_3=-\infty }^{\infty }\sum _{k_1}^{K_1-1}\sum _{k_2}^{K_2-1}\sum
   _{k_3}^{K_3-1}Q_j^{\mathit{ABC}}(\mathbf k,\mathbf n)\sum _{\mathbf m\neq 0}\left[\prod _{\mu
   }b_{\mu }(\mathbf m)e^{2\pi i\frac{m_{\mu }k_{\mu }}{K_{\mu }}}\right]f(\mathbf
   m)\\Q_j^{\mathit{ABC}}(\mathbf k,\mathbf n)=q_j\left[\frac{K_{\alpha }}{2\pi i}\right]^A\frac{\partial
   ^A}{\partial u_{\alpha j}^A}M_n(u_{\alpha j}-k_{\alpha }-n_{\alpha }K_{\alpha })\times \\\left[\frac{K_{\beta }}{2\pi
   i}\right]^B\frac{\partial ^B}{\partial u_{\beta j}^B}M_n(u_{\beta j}-k_{\beta }-n_{\beta }K_{\beta })\times
   \\\left[\frac{K_{\gamma }}{2\pi i}\right]^C\frac{\partial ^C}{\partial u_{\gamma j}^C}M_n(u_{\gamma j}-k_{\gamma
   }-n_{\gamma }K_{\gamma })\end{gathered}

where :math:`\omega` is the per-particle contribution for particle
:math:`j`, :math:`q` is the charge, other values are defined in the SPME
section (see: :ref:`SPME`)

File
~~~~

The heat flux method creates a file called HEATFLUX which contains the
relevant data structured as 6 Reals per line: STEP TEMPERATURE VOLUME HEAT-FLUX (x, y, z).

If **momentum_density** is specified with some atom types listed then the file 
will include additional data for each atom type's momentum density on each line 
after the HEAT-FLUX data. Given **momentum_density [Li]** in a Li-F simulation 
the data will be structured as 9 Reals per line and a Character: STEP TEMPERATURE VOLUME HEAT-FLUX (x, y, z) Li MOMENTUM-DENSITY-Li (x, y, z).


