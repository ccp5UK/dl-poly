.. _external-field:

External Fields
===============

.. index:: single: potential;external fields

In addition to the molecular force field, DL_POLY_4 allows the use of an
*external* force field. Examples of fields available include:

#. Electric field: (\ **elec**)

   .. math:: \underline{F_{i}} = \underline{F_{i}} + q_{i} \cdot \underline{E}

#. Oscillating shear: (\ **oshr**)

   .. math:: \underline{F}_{x} =A \cos(2n \pi \cdot z/L_{z})

#. Continuous shear: (\ **shrx**)

   .. math:: \underline{v}_{x} = \frac{1}{2} A \frac{|z|}{z}~~~~~:|z|>z_{0}

#. Gravitational field: (\ **grav**)

   .. math:: \underline{F_{i}} = \underline{F_{i}} + m_{i} \cdot \underline{G}

#. Magnetic field: (\ **magn**)

   .. math:: \underline{F_{i}} = \underline{F_{i}} + q_{i} \cdot (\underline{v_{i}} \times \underline{H})

#. Containing sphere: (\ **sphr**)

   .. math:: \underline{F}=A~(R_{0}-r)^{-n}~~~~~: r>R_{\rm cut}

#. Repulsive wall: (\ **zbnd**)

   .. math:: \underline{F}_{z}=A~(z_{o}-z)~~~~~: f \cdot z > f \cdot z_{o}~~,

   where :math:`f=+/-1` with default of :math:`1`.

#. X-Piston: (\ **xpis**)

   .. math::

      \underline{F}_{x} = \frac{m_{k}}{\sum_{k=i}^{j}m_{k}}
      P \cdot \underline{\textrm{Area}}(\perp\textrm{X-}direction)~~~~~: \forall~k=i,..,j~~.

#. Harmonic restraint zone in z-direction: (\ **zres**)

   .. math::

      \underline{F}_{z} = \left\{ \begin{array} {l@{\quad:\quad}l}
      A~(z_{com}-z_{max}) & z_{com} > z_{max} \\
      A~(z_{min}-z_{com}) & z_{com} < z_{min}~~,
      \end{array} \right.

   where :math:`z_{com}` is the chosen molecule centre of mass.

#. Harmonic restraint zone in z-direction :math:`-` (push out): (\
   **zrs\ :math:`-`**)

   .. math::

      \underline{F}_{z} = \left\{ \begin{array} {l@{\quad:\quad}l}
      A~(z-z_{max}) & z \ge (z_{max}+z_{min})/2 \\
      A~(z_{min}-z) & z < (z_{max}+z_{min})/2
      \end{array} \right.

#. Harmonic restraint zone in z-direction + (pull in): (\ **zrs+**)

   .. math::

      \underline{F}_{z} = \left\{ \begin{array} {l@{\quad:\quad}l}
      A~(z-z_{max}) & z > z_{max} \\
      A~(z_{min}-z) & z < z_{min}
      \end{array} \right.

#. Oscillating electric field: (\ **osel**)

   .. math:: \underline{F_{i}} = \underline{F_{i}} + q_{i} \cdot \underline{E} \cdot \sin(2 \pi \omega t)~,

   where :math:`t` is the simulated time.

#. Umbrella sampling (harmonic restraint)
   :cite:`torrie-77a,kastner-11a`: (\ **ushr**)

   .. math:: U_{AB} = \frac{k}{2} (R_{AB}-R_{0})^{2}~,

   is an umbrella sampling harmonic restraint between the centres of
   masses of two molecules (non-overlapping clusters of particles),
   :math:`A` and :math:`B`, with a force constant, :math:`k`, and an
   equilibrium distance, :math:`R_{0}`.

It is recommended that the use of an external field should be
accompanied by a :index:`thermostat` (this does not apply to examples 6 and 7,
since these are conservative fields). The “Oscillating shear” and
“X-piston” fields may only be used with orthorhombic cell geometry
(``imcon``\ :math:`=1`,:math:`2``) and “Continuous shear” field with slab cell geometry
(``imcon``\ :math:`=6`).

In the case of the “X-piston” field it is strongly advised that the
number of piston particles is chosen to be very small in comparison with
the rest of the system (:math:`< 5\%`) and that the piston contains its
own set of whole molecules (i.e. there are no molecules partially mapped
on the piston), which do not include any core-shell, CB, PMF or RB
units! The field releases the system’s centre of mass to move
unconstrained and gain momentum. This makes any temperature control
options control the full kinetic energy of the system and thus the only
ensemble valid under this conditions and possible within DL_POLY_4 at
the present is the micro-canonical (NVE)!

The user is advised to be careful with the parameters’ units! For more
insight, do examine
Table :numref:`(%s)<external-field-table>` and the
example at
equation :eq:`external-field-units_eq` in
Section :ref:`field-file`.

.. index:: single:potential;external field
   
In DL_POLY_4 external field forces are handled by the routines
``external_field_apply`` and .