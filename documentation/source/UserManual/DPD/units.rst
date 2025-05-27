.. _DPD_units:

DPD units
=========

DL_POLY_5 includes a unit scheme specifically for Dissipative Particle Dynamics (DPD) calculations,
that can be specified in both :ref:`CONTROL <control-file>` and :ref:`FIELD <field-file>` files. 
A directive in each of these files needs to be specified:

* The **io_units_scheme** directive in CONTROL needs to be set to **dpd**, so DL_POLY_5 will correctly
  output its results in DPD units. **Note that this option is not available for older-style CONTROL files!**
* The **units** directive in FIELD needs to be set to **dpd** so the user can provide interaction parameters 
  and particle properties (masses etc.) in DPD units.

Unlike the other unit schemes that can and have to be related to each other in *molecular* terms, 
the DPD units do not *need* to relate *directly* to atoms or molecules. This scheme relies instead
upon the connections between the different units and expresses all quantities in a consistent manner.
This is the same approach used by DL_MESO's DPD code, DL_MESO_DPD :cite:`seaton-13a`.

The three fundamental units in DL_POLY_5's DPD unit scheme are:

* Mass, :math:`[M]` (**dpd_m**)
* Length, :math:`[L]` (**dpd_l**)
* Energy, :math:`[E]` (**dpd_e**)

all of which can be equated to 'real' quantities for a given DPD simulation. Respectively, these are 
often defined as the mass of a specific bead :math:`m`, the cutoff distance for a given interaction 
:math:`r_{c}` and the thermal energy scale for the required system temperature :math:`k_{B} T`. 
Within DL_POLY_5's DPD units scheme, these quantities are *fictitiously but consistently* set to 
DL_POLY's internal units, i.e.

* :math:`[M] = 1` Dalton
* :math:`[L] = 1` Angstrom
* :math:`[E] = 10` J mol\ :sup:`-1`

Most of the other quantities can be expressed in terms of these three fundamental units. These and
their equivalents in DL_POLY internal units are shown in the table below.

.. table::
    :name: DPD Units Table

    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Quantity                |DPD unit name  |Equivalent in fundamental DPD units                 |Equivalent in DL_POLY internal units                                           |
    +========================+===============+====================================================+===============================================================================+
    |Mass :math:`[M]`        | **dpd_m**     | :math:`[M]`                                        | **internal_m**, :math:`m_o` = 1 Da                                            |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Length :math:`[L]`      | **dpd_l**     | :math:`[L]`                                        | **internal_l**, :math:`{\ell}_o` = 1 Å                                        |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Energy :math:`[E]`      | **dpd_e**     | :math:`[E]`                                        | **internal_e**, :math:`E_o` = 10 J mol\ :sup:`-1`                             |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Time :math:`[t]`        | **dpd_t**     | :math:`[L] \sqrt{\frac{[M]}{[E]}}`                 | **internal_t**, :math:`t_o` = 1 ps                                            |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Temperature :math:`[T]` | **dpd_temp**  | :math:`[E] / k_B`                                  | :math:`1/k_{B}` K (:math:`\sim 1.2027` K)                                     |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Pressure :math:`[P]`    | **dpd_p**     | :math:`[E] [L]^{-3}`                               | **internal_p**, :math:`{\cal P}_{o} \approx 16.61` MPa (163.9 atm)            |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Force :math:`[F]`       | **dpd_f**     | :math:`[E] [L]^{-1} \equiv [M] [L] [t]^{-2}`       | **internal_f** = 1 Da Å ps\ :sup:`-2`                                         |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Velocity :math:`[V]`    | **dpd_v**     | :math:`\sqrt{\frac{[M]}{[E]}} \equiv [L] [t]^{-1}` | **internal_v** = 1 Å ps\ :sup:`-1`                                            |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+
    |Charge :math:`[Q]`      | **dpd_q**     | :math:`[q]`                                        | **internal_q**, :math:`q_{o} = |e| \approx 1.602 \times 10^{-19}` C           |
    +------------------------+---------------+----------------------------------------------------+-------------------------------------------------------------------------------+

It should be noted that:

#.  The unit names in the above table should be used in the CONTROL file when specifying simulation parameters.
    For instance, cutoff distances, temperature and timestep can be expressed as:

    ::

      vdw_cutoff   1.0 dpd_l
      cutoff       3.0 dpd_l
      temperature  1.0 dpd_temp
      timestep     0.01 dpd_t
      
    With the exception of temperature, if any of these quantities are specified in DL_POLY's internal units
    directly (e.g. cutoffs in angstroms), they will still be assumed to be in DPD units.

#.  If **units dpd** is used in the second line of the FIELD file, all quantities in this file and other
    related files (including CONFIG and TABVDW files) will be read as though they are in DPD units. 
    That means all energy-based parameters will be in units of :math:`[E]`, all distances and positions 
    in :math:`[L]`, all masses in :math:`[M]` etc. This also extends to external field parameters, e.g.
    when applying gravitational, electric or magnetic fields. Since these units are implicitly used by 
    DL_MESO_DPD :cite:`seaton-13a`, FIELD and CONFIG files from this code can therefore be used by 
    DL_POLY_5 with little to no modification.

#.  The time unit :math:`[t] = [L] \sqrt{\frac{[M]}{[E]}}` is based upon matching up dimensions and units,
    and is sometimes referred to as the *natural* time unit. The alternative of matching the time unit to 
    a required diffusivity :cite:`Groot2001` (see below) can increase this value somewhat, but unlike the
    natural time unit it cannot be used directly with DL_POLY's molecular unit systems as its size can cause
    calculational instabilities and it would break other unit consistencies, e.g. for forces and velocities.

#.  The temperature used to define the energy scale :math:`[E] = k_{B} T` *often* equals the specified
    temperature for a DPD calculation and is set to 1, but this does not have to be the case and therefore 
    DL_POLY_5 can define the energy scale and temperature separately. Internally the DPD temperature unit 
    is equivalent to :math:`1/k_{B}` K to deal with several places in DL_POLY_5 where the Boltzmann constant 
    :math:`k_{B}` is used. 
    
#.  The DPD unit for charge :math:`[q]` is equivalent to DL_POLY_5's internal charge unit,  
    the absolute charge on an electron. Since the Coulombic conversion factor :math:`\gamma_{o}` scales 
    with both length and energy (:math:`{\ell}_{o}` and :math:`E_o`), we do not recommend directly 
    specifying the dielectric constant :math:`\epsilon` in the CONTROL file for DPD simulations. 
    Instead, we suggest specifying the Bjerrum length :math:`\lambda_B`, defined as

    .. math:: \lambda_B = \frac{q_{o}^2}{4 \pi k_{B} T \epsilon_{0} \epsilon}

    which can be scaled by :math:`[L]` for DPD simulations and set using the **coul_bjerrum_length**
    directive in the CONTROL file. The value :math:`\frac{4 \pi \lambda_B}{r_{c}}` is equivalent to 
    :math:`\Gamma`, another permittivity parameter often used for DPD simulations with 
    electrostatics :cite:`Groot2003`.

#.  If the CONTROL file specifies DPD units (**io_unit_scheme dpd**), all output files, including trajectories
    in HISTORY and system-wide properties given in OUTPUT and STATIS files, will be written in DPD units: 
    e.g. positions in :math:`[L]`, velocities in :math:`[V]`, forces in :math:`[F]`. 
    
    - Temperatures will be reported in terms of :math:`k_B T` or as though :math:`k_B = 1`, which is the 
      customary way temperature is defined in DPD simulations.
    
    - Pressures will be reported in the appropriate consistent units (:math:`[E] [L]^{-3}`) rather than the 
      default conversion to kilo-atmospheres [1]_.  

    - Interfacial tensions will also be reported in the appropriate consistent units (:math:`[F] [L]^{-1}`) 
      rather than converted to dynes per centimetre (:math:`10^{-3}` N m\ :sup:`-1`).
    
    - If this directive is not included, DL_POLY will (erroneously) report temperatures in Kelvin, pressures in
      kilo-atmospheres and interfacial tensions in dyn/cm. While this means the reported values will be incorrect,
      the accuracy of the DPD calculation itself will not be affected.
      
All of the quantities in the above table may alternatively be equated to 'real-life' 
values for a given DPD simulation. For instance, if a bead of water at room temperature (298.15K) 
consists of 3 molecules and the bead density for a DPD simulation is customarily set to 
:math:`3 [L]^{-3}`, the mass, length and energy units for this DPD simulation would be equivalent 
to:

* :math:`[M] = 54.04584` Da
* :math:`[L] = 6.46378` Å
* :math:`[E] = 2478.957` J mol\ :sup:`-1`

and it would be possible to run this simulation in DL_POLY_5 with everything (including derived 
units) scaled in the usual molecular units. The exception to this would be setting the DPD time 
unit :math:`[t]` to a value matching up the self-diffusivity obtained from simulations (e.g. 
measured using mean squared displacements, MSDs) with a real-life value :cite:`Groot2001`. 
For instance, if a water bead consists of :math:`N_w` molecules, the self-diffusivity of water
would equal:

.. math:: D_w = N_w D_{sim} [L]^2 [t]^{-1}

where :math:`D_{sim}` is the value obtained from simulations. For an example with :math:`N_w = 3`,
the natural time unit would be:

* :math:`[t] = 3.01810` ps

while the value obtained from matching diffusivities [2]_ would be:

* :math:`[t] \approx 93.1` ps

While the natural time unit could be used directly in DL_POLY_5's molecular units, the 
diffusivity-matched value cannot as the resulting timestep (a few picoseconds) may be too large 
to integrate particle forces without causing numerical instabilities. In this case, the 
simulation can only be carried out using DPD units.

.. [1]
   It should be noted that obtaining absolutely correct thermodynamic behaviour with typical DPD-based 
   interactions (e.g. Groot-Warren interactions) is very difficult, particularly since they are more 
   often parameterised for isothermal compressibility rather than pressure. As such, pressures in DPD 
   are normally only considered indicative for e.g. equilibration. They can be used to specify barostats 
   for e.g. NPT ensembles, although the required pressure for a DPD system might differ somewhat from the 
   value needed for an equivalent atomistic MD simulation!

.. [2]
   This follows the details given by Groot and Rabone :cite:`Groot2001`, which sets :math:`A_{ij} = 78 [E][L]^{-1}`,
   :math:`r_c = [L]`, and :math:`\gamma_{ij} = 4.5 [E][t][L]^{-2}`: the measured diffusivity
   was :math:`D_{sim} \approx 0.1707 [L]^{2} [t]^{-1}` and an experimental value for the self-diffusivity 
   of water at room temperature is :math:`D_w = 2.299 \times 10^{-9}` m\ :sup:`2` s\ :sup:`-1`.

