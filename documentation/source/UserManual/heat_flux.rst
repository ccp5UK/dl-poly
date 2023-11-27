.. _heat-flux:

Heat Flux
=========

Introduction
~~~~~~~~~~~~

It is possible to use DL_POLY_4 to calculate the heat flux of a material
with two-body interactions in an MD simulation. The heat flux can
subsequently be used as a means to calculate the thermal conductivity of
a material via the Green-Kubo relation. Currently, the heat flux is only
viable for two-body interactions, but valid for any two-body
interactions.

To enable the calculation of heat flux add the **heat_flux** keyword
into the CONTROL file.

Theory
~~~~~~

The heat flux for two-body interactions is defined as:

.. math:: \underline{J} = \frac{1}{V} \left[ \sum\limits^{N}_{i} e_{i} \underline{v}_{i} - \sum\limits^{N}_{i} \underline{\underline{\textbf{S}}}_{i} \underline{v}_{i} \right]\label{heat-flux-definition}
   :label: heat-flux-definition_eq

where :math:`\underline{\underline{\textbf{J}}}` is the heat flux,
:math:`V` is the volume of the cell, :math:`N` is the number of
particles, :math:`e` is the energy, :math:`\underline{v}` is the
velocity, :math:`\underline{\underline{\textbf{S}}}` is the stress. All
subscript :math:`i` refer to the particle.

The thermal conductivity can then be as an auto-correlation of the heat
flux over a run:

.. math:: \kappa = \frac{V}{k_{B} T^{2}} \int\limits_{0}^{\infty} \langle \underline{J}(0)  \underline{J}(t) \rangle \, \mathrm{d}t

where :math:`\kappa` is the thermal conductivity, :math:`V` is the
volume of the cell, :math:`k_{B}` is the Boltzmann constant, :math:`J`
is the heatflux.

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
relevant data structured as: STEP PRESSURE VOLUME HEAT-FLUX
