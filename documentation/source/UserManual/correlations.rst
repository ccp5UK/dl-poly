.. _correlation-functions:

Correlation Functions
=====================

Introduction
~~~~~~~~~~~~

DL_POLY_5 includes on the fly computation of time correlations of
observable quantities. This functionality allows for key correlation
functions, and derived quantities, to be calculated without saving
trajectory data.

The framework has been designed to support arbitrary correlation
function. Currently, for the user, implemented observable quantities
include: Atom velocity and the system’s stress and heat flux (with
two-body interactions, refer to Section :ref:`heat-flux`).
Any pair of these observables can be composed into a correlation
function. However a typical use case is to compute the velocity (VAF),
stress (SAF), and heat flux auto-correlations (HFAF). These allow the
computation of the vibrational density of states (via a Fourier
transform of the VAF) or shear-viscosity and thermal conductivity from
Green-Kubo relations on the latter two. In fact the shear-viscosity and
thermal conductivity are calculated whenever the required correlation
functions are requested by the user in the (new style) CONTROL file (see
Section :ref:`new-control-file`). For detailed usage instructions
refer to the User Control in Section :ref:`cor-user`.

Theory
~~~~~~

Taking as an example the shear-stress auto-correlation function (in
discrete form),

.. math:: C^{\tau} = \frac{1}{T-\tau}\sum_{t'}^{T} \sigma_{xy}^{t'}\sigma_{xy}^{t'+\tau} 
   :label: stress-cor_eq

where :math:`\tau` indicates a discrete lag time, and
:math:`\sigma_{xy}^{t}` indicates the :math:`xy` component of stress at
discrete simulation time step :math:`t`. Using a STATIS file (see
Section :ref:`statis-file`) the values of
:math:`\sigma_{xy}^{t}` can be read post simulation and calculated
directly. This method however requires either saving stress data for all
discrete time points :math:`1,2, \ldots T`, or suffering a loss of
accuracy by sub-sampling them. Additionally with correlation functions
such as the velocity auto-correlation function, per-atom data is also
required. Long timescale simulations and/or large system sizes present a
scaling issue for both memory and run-time.

DL_POLY_5 utilises the Multiple-tau correlator
:cite:`Ramirez2010Efficient`, which is one particular on the
fly correlation algorithm that addresses these issues. Briefly the
method works by accumulating data in a series of hierarchical block
averages. Three parameters control this **correlation_blocks**,
**correlation_block_points**and **correlation_window** (written in
terms of CONTROL directives). These control the number of hierarchical
blocks, the number of distinct points within each block, and the length
of an averaging window between blocks respectively. A fourth parameter 
**correlation_update_frequency** controls how frequently each correlation 
is updated (the default is to follow **stats_frequency**).

In more detail, given an empty correlator, as new data is submitted to
it a sum is accumulated and the data points held in temporary storage.
When the first block remains unfilled (less than
**correlation_block_points** data points have been submitted) the
product of the new data point with each temporarily stored value is
added to a correlation accumulator for this first block. At the first
block all multiplications must be carried out, in subsequent blocks only
points between
:math:`\textbf{correlation\_block\_points}/\textbf{correlation\_window}`
and **correlation_block_points** need be updated. Once the first block
contains **correlation_block_points** data entries the sum divided by
the **correlation_window** is passed to the next level, and the
temporary data stored in the first level is cleared along with its
accumulated sum (but not its accumulated correlation). In this way the
complete correlation is accumulated over a system’s trajectory, with
only ephemeral storage of “raw” data. As the simulation progresses, past
values are retained at ever decreasing resolution replaced by current
data. In terms of storage complexity the algorithm scales as
:math:`(p-1)w^b` per correlation, where for brevity :math:`p` is
**correlation_block_points**:math:`w` is **correlation_window** and
:math:`b` is **correlation_blocks**. It should be noted that certain
correlations, such as velocity, are computed on a per-particle basis.
This requires :math:`N` correlators for a system of :math:`N` atoms.

Particular correlation functions can be used to analyse the results of a
simulation and compare to experimental systems. For example the SAF can
be integrated to yield a Green-Kubo relation for sheer-viscosity. That
is with analytic expressions for sheer-stress :math:`\sigma_{xy}(t)` in
continuous time :math:`t`, sheer-viscosity is

.. math:: \eta = \frac{V}{k_{b}T}\int_{0}^{\infty}dt' \langle \sigma_{xy}(0)\sigma_{xy}(t')\rangle.
   :label: viscosity-gk_eq

Where :math:`V` and :math:`T` are the system volume and temperature
respectively, with :math:`k_{b}` Boltzmann’s constant. The integral can
be performed numerically over the discretised form of the correlation
function in Equation :eq:`stress-cor_eq` to estimate
sheer-viscosity from simulation data. Similar relations exist for e.g.
HFAF and thermal conductivity i.e.

.. math:: \lambda = \frac{V}{3k_{b}T^2}\int_{0}^{\infty} dt' \langle \textbf{J}(0)\cdot \textbf{J}(t') \rangle. 
   :label: thermal-conductivity-gk_eq

The prefactor includes a multiplication with volume due to the
definition of heat flux, :math:`\textbf{J}(t)`, in DL_POLY_4 already
including a volume division, see
Equation :eq:`heat-flux-definition_eq`.

.. _cor-user:

User Control
~~~~~~~~~~~~

Input
^^^^^

In the (new style) CONTROL (see Section :ref:`new-control-file`)
correlations are specified by an array of observable pairs in the format
**A_CA-B_CB** where **A** and **B** are observables and **CA** and **CB** 
are components of those observables. These values may be STATIS records 1-27, e.g. ``volume``, (see
Section :ref:`statis-records`). Or may take the forms listed in 
Table :numref:`(%s)<tab-cor-control>` where long and short form can be mixed. 
Correlation options can be specified for each pair of observables separately, these
options are listed in Table :numref:`(%s)<tab-cor-options>`

.. _tab-cor-options:

.. table:: 
      User control directives for correlation options. All are vector options.

   ================================ ============== ===============================================
   Option                           Type           Purpose                  
   ================================ ============== ===============================================
   **correlation_observable**       String Vector  Set correlation observables See Table :numref:`(%s)<tab-cor-control>`
   **correlation_block_points**     Integer Vector Set correlation points per block
   **correlation_blocks**           Integer Vector Set correlation blocks
   **correlation_window**           Integer Vector Set averaging between block
   **correlation_update_frequency** Integer Vector Set update frequency in steps
   ================================ ============== ===============================================

For example to compute the VAF 
(across all dimensions) and SAF along just **xy** one may write,

::

       correlation_observable [v_x-velocity_x v_y-v_y v_z-v_z s_xy-stress_xy]
       correlation_block_points [600 600 600 5000]
       correlation_blocks [2 2 2 1]
       correlation_window [2 2 2 1]

DL_POLY_5 will then accumulate the VAF (x,y,z) with 2 blocks each with 600
points per block, and a averaging length 2 and the SAF (just xy) with a single
block with 5000 points. By default **correlation_window** :math:`=1`,
**correlation_block_points** :math:`=100` and **correlation_blocks**
:math:`=1`.

.. _tab-cor-control:

.. table:: 
      User control directives in the new style CONTROL file for on
      the fly correlations. Any combination of observables can be correlated.
      Except for **lc**, **tc**, **kd**, and **ec** which may only be correlated
      between each other. Observables indicated as Per-particle require computation, 
      and storage of data scaling with system size :math:`N`, as one correlator for 
      each atom is created. Correlators marked as Per-species are split automatically
      across atomic species. For current correlations the correlated :math:`{\bf k}`-space
      currents are averaged over atoms among each species before correlation. Observables 
      marked with "scalar" in components require no _x component specification.

   ========================= ========== ========================================== ============ ============ ======================================= =========================================
   String                    Short-hand Observable                                 Per-particle Per-species  Components                              Notes
   ========================= ========== ========================================== ============ ============ ======================================= =========================================
   **velocity**              **v**      particle velocity                          yes          yes          **x, y, z**
   **stress**                **s**      system stress                              no           no           **xx, xy, xz, yx, yy, yz, zx, zy, zz**
   **heat_flux**             **hf**     system heat flux                           no           no           **x, y, z**
   **longitudinal_current**  **lc**     Longitudinal component of current          no           yes          **x, y, z**                             Automatically calculated over KPOINTS
   **transverse_current**    **tc**     Transverse component of current            no           yes          **x, y, z**                             .
   **kdensity**              **kd**     k-space density                            no           yes          scalar                                  .
   **energy_density**        **ed**     energy current                             no           yes          **x, y, z**                             .
   **energy_current**        **ec**     energy current                             no           yes          **x, y, z**                             ., requires **energy_stress_currents On**
   **kstress**               **ks**     energy current                             no           yes          **xx, xy, xz, yx, yy, yz, zx, zy, zz**  ., requires **energy_stress_currents On**
   ========================= ========== ========================================== ============ ============ ======================================= =========================================

The maximum lag time for a correlation is given by 

.. math:: (p-1)m^{l}\Delta t,
      :label: cor_max_lag

where :math:`p`` is **correlation_block_points** and :math:`m` 
is **correlation_window**and :math:`l` is **correlation_blocks**. 
The algorithmic scaling is bounded above by,

.. math:: p\frac{m+1}{m}.
      :label: cor_scaling

For per-atom correlations such as velocity updating one 
correlator will also scale with an additional factor of :math:`N`, 
the atom count.

These facts can be combined to obtain a given correlation timescale
whilst balancing performance. I.e. an higher :math:`p` increase accuracy,
but has the greatest performance impact. Increase :math:`m` and :math:`l`
have comparatively less performance impact but give access to longer 
correlation lengths at some reduced accuracy.

Correlating STATIS Values
^^^^^^^^^^^^^^^^^^^^^^^^^

All values in the default STATIS output (records 1-27, see Section :ref:`statis-records`)
may be correlated with any other correlation observable. These are listed in Table 
:numref:`(%s)<tab-statis-cor>`. Each are scalar values with no component specification 
needed.

.. _tab-statis-cor:

.. table::
      Observable STATIS valued which can be correlated.

   ============ ==========
   String       Observable
   ============ ==========
   **eng_tot**  total extended system energy, :math:`E^{x}_{tot}=(E_{kin}+E_{rot})+E_{conf}+E_{consv}` (i.e. including the conserved quantity, :math:`E_{consv}`)
   **temp_tot** system temperature, :math:`2\frac{E_{kin}+E_{rot}}{f k_{B}}`
   **eng_cfg**  configurational energy, :math:`E_{conf}`
   **eng_src**  short range potential energy
   **eng_cou**  electrostatic energy
   **eng_bnd**  chemical bond energy
   **eng_ang**  valence angle and 3-body potential energy
   **eng_dih**  dihedral, inversion, and 4-body potential energy
   **eng_tet**  tethering energy
   **eng_pv**   enthalpy (:math:`E^{x}_{tot} + {\cal P} \cdot V`) for NVE/T/E\ :math:`_{kin}` ensembles enthalpy (:math:`E^{x}_{tot} + P \cdot {\cal V}`) for NP/\ :math:`\sigma`\ T or NP\ :math:`_{n}`\ A/\ :math:`\gamma` ensembles
   **temp_rot** rotational temperature, :math:`E_{rot}`
   **vir_cfg**  total virial
   **vir_src**  short-range virial
   **vir_cou**  electrostatic virial
   **vir_bnd**  bond virial
   **vir_ang**  valence angle and 3-body virial
   **vir_con**  constraint bond virial
   **vir_tet**  tethering virial
   **volume**   volume, :math:`{\cal V}`
   **temp_shl** core-shell temperature
   **eng_shl**  core-shell potential energy
   **vir_shl**  core-shell virial
   **alpha**    MD cell angle :math:`\alpha`
   **beta**     MD cell angle :math:`\beta`
   **gamma**    MD cell angle :math:`\gamma`
   **vir_pmf**  PMF constraint virial
   **press**    pressure, :math:`{\cal P}`
   ============ ==========

Currents
^^^^^^^^

:math:`{\bf k}`-space currents (and density) may be calculated by specifying a KPOINTS file 
see Sections :ref:`currents` and :ref:`kpoints-file_sec`. As shown in 
Table :numref:`(%s)<tab-cor-control>` a user may auto- or cross-correlate all of these values 
between eachother. Like other correlations these observables must be requested component 
wise, e.g. :math:`x, y, z`. However, DL_POLY_5 will automatically generate correlators for all
user specified KPOINTS across all distinct atom types. 

For example with a Lithium-Flouride (LiF) 
system and 3 user KPOINTS requesting the correlation **lc_x-tc_z** will calculate the 
cross correlations of the :math:`x`` component of the longitudinal momentum current 
with the :math:`z` component of the transverse momentum current for all 3 KPOINTS for both
Li and F separately. This will generate 6 correlations in COR.

Current correlations will have output names of the form **atomtype-shortname_kpointnumber_component** 
for example **Li-ks_1_zz-Li-ks_1_zz** and **Li-kd_10-Li-kd_10**. See also the output section below.

Complex Observables
^^^^^^^^^^^^^^^^^^^

Any complex observable is correlated according to the formula 

.. math:: C(A, B) = \langle A\cdot \text{Conj}(B) \rangle,

where :math:`\text{Conj}(x)` is the complex conjugate of :math:`x`.

Output
^^^^^^

When correlation functions are specified by the user the resulting data
is written as a YAML file, COR, containing the distinct correlations
with their lag times, components. The header section may contain derived 
data computed from the correlation functions available. For
example when computing the SAF the viscosity and kinematic viscosity 
value for the simulation is automatically calculated. The value will be
an averaged over commensurate correlations, i.e. if **s_xy-s_xy** and **s_yz-s_yz** 
are specified they will be both be used to compute and average viscosity. The 
individual components comprising the observable data (if applicable) will be 
outputted in a **components** array of the observables section.
An example COR file is:

::

      %YAML 1.2
      ---
      title: 'CONFIG generated by ASE'
      observables:
            viscosity:
                  value:   0.24150069E-03
                  components: [  0.39309725E-03,  0.89904136E-04]
                  units: Katm ps 
            kinematic-viscosity:
                  value:   0.21993727E-03
                  components: [  0.35799788E-03,  0.81876662E-04]
                  units: Katm ps / (amu / Ang^3)
            thermal-conductivity:
                  value:   0.96869183E-06
                  units: e.V / (ps Ang K)
            elasticity_tensor:
                components: [C_xxxx , C_xxyy , C_xxzz , C_yyyy , C_yyzz , C_zzzz , C_yzyz , C_zxzx , C_xyxy]
                values: [   24.007901    ,   10.219979    ,   11.545088    ,   17.599676    ,   8.3646309    ,   19.861066    ,   14.510356    ,   16.062736    ,   14.391474    ]
                units: Katm
      correlations:
            stress_xy-stress_xy:
                  parameters:
                        points_per_block: 5000
                        number_of_blocks: 1
                        window_size: 1
                  lags: [   0.0000000    ,   1.0000000    ,   ... ]
                  value: [  0.28183438E-03,  0.28071854E-03,  ... ]
            stress_yz-stress_yz:
                  parameters:
                        points_per_block: 100
                        number_of_blocks: 1
                        window_size: 1
                  lags: [   0.0000000    ,   1.0000000    ,   ... ]
                  value: [  0.64655167E-04,  0.64642762E-04,  ... ]
            heat_flux_x-heat_flux_x:
                  parameters:
                        points_per_block: 100
                        number_of_blocks: 1
                        window_size: 1
                  lags: [   0.0000000    ,   1.0000000    ,   ... ]
                  value: [  0.23606349E-08,  0.23606349E-08,  ... ]
            Ar-velocity_x-velocity_y:
                  parameters:
                        points_per_block: 100
                        number_of_blocks: 1
                        window_size: 1
                  lags: [   0.0000000    ,   1.0000000    ,   ... ]
                  value: [  0.14340287E-01,  0.14342140E-01,  ... ]

Specific Correlation Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _tab-cor-observables:

.. table:: 
      Possible observable quantities written out to the COR file.

   ======================== ====================== =================
   Name                     Related correlation(s) Description        
   ======================== ====================== ================= 
   **viscosity**            shear-stress.          From Green-Kubo formula :eq:`viscosity-gk_eq`.
   **kinematic-viscosity**  shear-stress.          Viscosity divided by density.
   **thermal-conductivity** heatflux, diagonal.    From Green-Kubo formula :eq:`thermal-conductivity-gk_eq`.
   **elasticity-tensor**    stress.                From stress-fluctuation method.
   ======================== ====================== =================

**Stress correlations**\ : When accumulating system stress correlation functions a derived sheer-viscosity 
measurement is written to the output file, along with a kinematic viscosity. These are both calculated using equation 
:eq:`viscosity-gk_eq`, averaged over any of the xy, yz, zx, yx, zy, and xz  correlations requested. 
The latter is calculated by additionally dividing by the system density, thus the kinematic viscosity is 

.. math::

    \eta_k = \eta / \rho,

where :math:`\rho` is averaged over the simulation time.

**Heatflux correlations**\ : When accumulating heatflux correlations the thermal-conductivity is written to output as a 
derived measurement. This is calculated following equation :eq:`thermal-conductivity-gk_eq` (i.e. with averaging over x, y, 
and z directions if any correlations are specified).
