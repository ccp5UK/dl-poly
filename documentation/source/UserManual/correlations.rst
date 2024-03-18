Correlation Functions
=====================

Introduction
~~~~~~~~~~~~

DL_POLY_4 includes on the fly computation of time correlations of
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

DL_POLY_4 utilises the Multiple-tau correlator
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
**x-y** where **x** and **y** may take the string values in
Table :numref:`(%s)<tab-cor-control>`. For example to compute
the VAF and SAF one may write,

::

       correlation_observable [velocity-velocity s-s]
       correlation_block_points [600 5000]
       correlation_blocks [2 1]
       correlation_window [2 1]

DL_POLY_4 will then accumulate the VAF with 2 blocks each with 600
points per block, and a averaging length 2 and the SAF with a single
block with 5000 points. By default **correlation_window** :math:`=1`,
**correlation_block_points** :math:`=100` and **correlation_blocks**
:math:`=1`.

.. _tab-cor-control:

.. table:: 
      User control directives in the new style CONTROL file for on
      the fly correlations. Any combination of observables can be correlated.
      Observables indicated as per-particle require storage of data scaling
      with system size :math:`N`, as one correlator for each atom is created.

   ============= ========== ================= ============
   String        Short-hand Observable        Per-particle
   ============= ========== ================= ============
   **velocity**  **v**      particle velocity yes
   **stress**    **s**      system stress     no
   **heat_flux** **hf**     system heat flux  no
   ============= ========== ================= ============


Through a combination of the three parameters short or long timescale
correlations may be computed. For example taking **correlation_blocks**
:math:`= 1`, **correlation_block_points** :math:`= 100`, and
**correlation_window** :math:`= 1` will result in a maximum lag time
correlated of :math:`99 \Delta t` for simulation time step
:math:`\Delta t`. Whereas taking **correlation_window** :math:`= 2` and
**correlation_blocks** :math:`= 2` will give a higher maximum
correlation lag time, :math:`(198 \Delta t)`, but the averaging window
of :math:`2` will result in a small accuracy reduction.

Output
^^^^^^

When correlation functions are specified by the user the resulting data
is written as a YAML file, COR, containing the distinct correlations
with their lag times, components, and any derived quantities. For
example when computing the SAF the derived viscosity value for the
simulation is automatically calculated, in this case the COR output file
may look like the following

::

   %YAML 1.2
   ---
   title: argon fcc initial conditions
   correlations:
       - name: [stress-stress                    , global]
         parameters:
               points_per_block: 5000
               number_of_blocks: 1
               window_size: 1
         derived:
               viscosity:
                     value:   0.57490361    
                     units: Katm ps 
         lags: [   0.0000000    ,  0.10000000E-03,  0.20000000E-03, ...]
         components: 
              stress_xx-stress_xx: [   1.0484384,   1.0484383,   1.0484382, ...]
              stress_xy-stress_xy: [                   ...                     ]
              ...
              stress_zz-stress_zz: [                   ...                     ]


Specific Correlation Output 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Stress correlations**\ : When accumulating system stress correlation functions a derived sheer-viscosity 
measurement is written to the output file, along with a kinematic viscosity. These are both calculated using equation 
:eq:`viscosity-gk_eq`, averaged over the xy, yz, and zx correlations. The latter is calculated by additionally dividing 
by the system density, thus the kinematic viscosity is 

.. math::

    \eta_k = \eta / \rho,

where :math:`\rho` is averaged over the simulation time.

**Heatflux correlations**\ : When accumulating heatflux correlations the thermal-conductivity is written to output as a 
derived measurement. This is calculated following equation :eq:`thermal-conductivity-gk_eq` (i.e. with averaging over x, y, 
and z directions). Additionally the components of lattice thermal-conductivity are also written for each component pairs, 
xx, xy, xz, etc.
