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

.. math:: C^{\tau} = \frac{1}{T-\tau}\sum_{t'}^{T} \sigma_{xy}^{t'}\sigma_{xy}^{t'+\tau} \label{stress-cor}
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
**correlation_block_points**, and **correlation_window** (written in
terms of CONTROL directives). These control the number of hierarchical
blocks, the number of distinct points within each block, and the length
of an averaging window between blocks respectively.

In more detail, given an empty correlator, as new data is submitted to
it a sum is accumulated and the data points held in temporary storage.
When the first block remains unfilled (less than
**correlation_block_points** data points have been submitted) the
product of the new data point with each temporarily stored value is
added to a correlation accumulator for this first block. At the first
block all multiplications must be carried out, in subsequent blocks only
points between
:math:`\textbf{correlation_block_points}/\textbf{correlation_window}`
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
**correlation_block_points**, :math:`w` is **correlation_window** and
:math:`b` is **correlation_blocks**. It should be noted that certain
correlations, such as velocity, are computed on a per-particle basis.
This requires :math:`N` correlators for a system of :math:`N` atoms.

Particular correlation functions can be used to analyse the results of a
simulation and compare to experimental systems. For example the SAF can
be integrated to yield a Green-Kubo relation for sheer-viscosity. That
is with analytic expressions for sheer-stress :math:`\sigma_{xy}(t)` in
continuous time :math:`t`, sheer-viscosity is

.. math:: \eta = \frac{V}{k_{b}T}\int_{0}^{\infty}dt' \langle \sigma_{xy}(0)\sigma_{xy}(t')\rangle.\label{viscosity-gk}
      :label: viscosity-gk_eq

Where :math:`V` and :math:`T` are the system volume and temperature
respectively, with :math:`k_{b}` Boltzmann’s constant. The integral can
be performed numerically over the discretised form of the correlation
function in Equation :eq:`stress-cor_eq` to estimate
sheer-viscosity from simulation data. Similar relations exist for e.g.
HFAF and thermal conductivity i.e.

.. math:: \lambda = \frac{V}{3k_{b}T^2}\int_{0}^{\infty} dt' \langle \textbf{J}(0)\cdot \textbf{J}(t') \rangle. \label{thermal-conductivity-gk}
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
are components of those observables. These values may take the 
forms listed in Table :numref:`(%s)<tab-cor-control>` where long and 
short form can be mixed. For example to compute the VAF 
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
      Observables indicated as per-particle require storage of data scaling
      with system size :math:`N`, as one correlator for each atom is created.

   ============= ========== ================= ============ =======================================
   String        Short-hand Observable        Per-particle Components
   ============= ========== ================= ============ =======================================
   **velocity**  **v**      particle velocity yes           **x, y, z**
   **stress**    **s**      system stress     no            **xx, xy, xz, yx, yy, yz, zx, zy, zz**
   **heat_flux** **hf**     system heat flux  no            **x, y, z**
   ============= ========== ================= ============ =======================================

The maximum lag time for a correlation is given by 

.. math:: (p-1)m^{l}\Delta t,
      :label: cor_max_lag

where :math:`p`` is **correlation_block_points** and :math:`m` 
is **correlation_window**, and :math:`l` is **correlation_blocks**. 
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
