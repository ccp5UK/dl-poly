Exdenting DL_POLY_4 to reactive systems: the Empirical Valence Bond method
==========================================================================

Framework and motivation
------------------------

A key component of DL_POLY_4 is the Force Field (FF) to model the interactions between atoms. As already described in previous chapters, such atomic interactions are often modelled by relatively simple functional forms with parameters either fitted to experimental data or derived from quantum mechanical calculations. In most of the classical FFs available, functional forms and fitted parameters remain unchanged during the course of the molecular dynamics (MD) simulation. Indeed, this is the type of FFs that DL_POLY_4 can handle. In reactive processes, however, the nature of the interactions inevitable changes due to the formation of new chemical species. For this reason, standard FFs (thence DL_POLY_4) are not suitable to simulate chemical reactions. We shall refer to such FFs as non-reactive.

An alternative to simulate chemical reactions is offered by the so called Reactive FFs (RFFs). In contrast to standard FFs where interactions are modelled for a particular state with a given topology and chemistry, RFFs are designed to model the interatomic interactions valid for multiple states that are chemically different. The task of designing RFFs, however, is very challenging and requires a high level of expertise to tackle a multi-dimensional problem, where the modelled interactions are often expressed by complicated functional forms with many strongly coupled parameters that are optimised via the use of sophisticated tools. Even though RFFs have evolved considerably in the last years, a general parametrization is not yet available and, instead, parameters have to be tuned to specific chemical systems and environments.

Within this framework and to the purpose of extending the applicability of DL_POLY_4 to simulate reactive processes, the Empirical Valence Bond (EVB) method :cite:`duarte2017` offers an appealing alternative for computational implementation and development. In contrast to facing the challenges of building RFFs, the EVB method defines a suitable matrix using computed quantities of the participating chemical states, where each state is modelled by a non-reactive FF. Via the definition of appropriate coupling terms and matrix diagonalization at each time step, it is possible to obtain potential energy landscapes that account for the change in chemistry when sampling conformations between the participating, chemically different, states.

In contrast to RFFs, the advantage of the EVB method lies in the large availability of standard non-reactive FFs libraries. In addition, despite the initial task to calibrate the coupling terms against reference data, research has demonstrated that these couplings are invariant to the surrounding electrostatics, making it possible to simulate the same reactive unit in different environments. For further details about the applications of the EVB method, we refer the user to ref. :cite:`scivetti-evb`.

The fundamentals of the EVB method are presented in the next section. Strategies to calibrate EVB-FFs are discussed in section :ref:`evb-calibrate`. The computational implementation of the EVB method is described in section :ref:`implement`. Finally, section :ref:`evb-users` provides a guideline to users on how to prepare the settings for EVB simulations with DL_POLY_4.  

The EVB Method 
--------------
.. _evbMethod:

Let us assume an atomic system composed of :math:`N_{p}` particles with positions described by the set  of vectors :math:`\mathbf{R}`. The non-reactive force field (FF) for the chemical state $m$ is described by the configurational energy :math:`E_{c}^{(m)}(\mathbf{R})` and the set of forces :math:`\vec{F}_{J}^{(m)}(\mathbf{R})`, where the index :math:`J` runs over the total number of particles. The configurational energy function :math:`E_{c}^{(m)}(\mathbf{R})` has the decomposition of eq. :eq:`decomp-ene`. In the following, however, we shall omit the presence of external fields, such as electric or magnetic. In the current notation, we shall use indexes :math:`m` and :math:`k` for the chemical states (and FFs), :math:`I` and :math:`J` for atoms and Greek letters for Cartesian coordinates. Indexes in parenthesis are used to emphasize the particular chemical state.

The purpose of the EVB method is to couple :math:`N_F` non-reactive force fields to obtain a reactive potential. These FFs are coupled through the Hamiltonian :math:`\hat{H}_{\text{EVB}}` with a matrix representation :math:`H_{\text{EVB}} \in \mathcal{R}^{N_F \times N_F}` that has the following components

.. _evbmatrix:
.. math:: 
   :label: evbmatrix_eq

    H^{mk}_{\text{EVB}}(\{\bf{R}\})=\begin{cases} E_{c}^{(m)}(\{\mathbf{R}\})               \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,  m=k   \\
    C_{mk}(\epsilon_{mk})                     \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,   m \ne k 
    \end{cases}

where each diagonal element corresponds to the configurational energy :math:`E_{c}^{(m)}(\mathbf{R})` of the non-reactive FF that models the interactions as if the system was in the chemical state :math:`(m)`, whereas the off-diagonal terms C\ :math:`_{mk}` are the couplings between states :math:`m` and :math:`k`. For convenience in the notation, we shall omit hereinafter the dependence on the set of coordinates :math:`\mathbf{R}` for the particles. Even though there are different possible choices for the coupling terms, in the above definition we have set :math:`C_{mk}` to depend on :math:`\epsilon_{mk}=E_{c}^{(m)}-E_{c}^{(k)}=-[E_{c}^{(k)}-E_{c}^{(m)}]=-\epsilon_{km}`, where :math:`\epsilon_{mk}` is commonly referred to as energy gap and defines a possible reaction coordinate for the reactive process :cite:`mones2009`. Since the :math:`H_{\text{EVB}}` matrix is Hermitian by construction and the :math:`C_{mk}` terms are real, the condition of :math:`C_{mk}=C_{km}` must be imposed to the off-diagonal elements. Diagonalization of :math:`H_{\text{EVB}}` leads to :math:`N_F` possible eigenvalues :math:`\{\lambda_1,...,\lambda_{N_{F}}\}` with

.. math:: 
    
    H_{\text{EVB}}\Psi_{\lambda_m}=\lambda_m \Psi_{\lambda_m}, \,\,\,\,\,\,\,\,\, m=1,...,N_F.

The EVB energy, :math:`E_{\text{EVB}}`, is defined as the lowest eigenvalue

.. math::
   :label: Eevb_eq

   \label{eq:Eevb}
       E_{\text{EVB}}=min(\lambda_1,...,\lambda_{N_F})

with the corresponding normalized EVB eigenvector

.. math::
   :label: Psi-evb-norm_eq

   \label{eq:Psi-evb-norm}
       \Psi_{\text{EVB}}=\Psi_{min(\lambda_1,...,\lambda_{N_F})}.

and

.. math::
   :label: EevbPsi_eq

   \label{eq:EevbPsi}
       E_{\text{EVB}}=\big\langle \Psi_{\text{EVB}}\big|\hat{H}_{\text{EVB}}\big| \Psi_{\text{EVB}}\big \rangle.

Since the eigenvector :math:`\Psi_{\text{EVB}}` is real and normalized
we have

.. math::

   \label{eq:evbPsinorm}
      \sum_{k=1}^{N_F} \big|\Psi^{(k)}_{\text{EVB}}\big|^{2}=1

from which we can interpret :math:`|\Psi^{(k)}_{\text{EVB}}\big|^{2}` as
the fraction of the chemical state :math:`(k)` being part of the EVB
state. The eigenvector :math:`\Psi_{\text{EVB}}` can also be represented
as a column vector :math:`\in \mathcal{R}^{N_F \times 1}` where
:math:`\Psi^{(k)}_{\text{EVB}}` is the element of the :math:`k`-row.
Thus, eq. :eq:`EevbPsi_eq` is expressed as a matrix
multiplication

.. math::
   :label: EevbPsimat_eq

   \label{eq:EevbPsimat}
      E_{\text{EVB}}=\sum_{m,k=1}^{N_F} \tilde{\Psi}^{(m)}_{\text{EVB}} H^{mk}_{\text{EVB}}\Psi^{(k)}_{\text{EVB}}

where :math:`\tilde{\Psi}_{\text{EVB}}` is the transpose of
:math:`{\Psi}_{\text{EVB}}`. The resulting EVB force over the particle
:math:`J`, :math:`\vec{F}_{J}^{\text{EVB}}`, follows from the
Hellman-Feynman theorem

.. math::
   :label: Fevb_eq

   \begin{aligned}
   \label{eq:Fevb}
      &\vec{F}_{J}^{\text{EVB}}=-\nabla_{\vec{R}_J}E_{\text{EVB}}=-\big\langle \Psi_{\text{EVB}}\big| \nabla_{\vec{R}_J} \hat{H}_{\text{EVB}} \big| \Psi_{\text{EVB}}\big \rangle \nonumber \\
      &= \sum_{\alpha=x,yz} F_{J\alpha}^{\text{EVB}} \,\, \check{\alpha}
      \end{aligned}

where :math:`\check{\alpha}` corresponds to each of the orthonormal
Cartesian vectors and

.. math::
   :label: Fevb2_eq

   \label{eq:Fevb2}
      F_{J\alpha}^{\text{EVB}}=-\big\langle \Psi_{\text{EVB}}\big| \frac{\partial \hat{H}_{\text{EVB}}}{\partial_{R_{J\alpha}}}\big| \Psi_{\text{EVB}}\big \rangle.

From eq. :eq:`evbmatrix_eq` the matrix components of
the operator
:math:`\frac{\partial \hat{H}_{\text{EVB}}}{\partial_{R_{J\alpha}}}` are
given as follows

.. math::
   :label: gradevb_eq

   \label{eq:gradevb}
      \frac{\partial H^{mk}_{\text{EVB}}}{\partial R_{J\alpha}}
      =\begin{cases}
      \frac{\partial E_{c}^{(m)}}{\partial R_{J\alpha}}=-F^{(m)}_{J\alpha} \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, m=k  \\
      \\
      \begin{aligned}
      \frac{d C_{mk}}{\partial R_{J\alpha}} &=\frac{d C_{mk}(\epsilon_{mk})}{d\epsilon_{mk}}\frac{\partial \epsilon_{mk}}{\partial R_{J\alpha}}\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,  m \ne k\\
                                                                &=\frac{d C_{mk}(\epsilon_{mk})}{d\epsilon_{mk}} \left[\frac{\partial E_{c}^{(m)}}{\partial J\alpha}-\frac{\partial E_{c}^{(k)}}{\partial J\alpha}\right]\\
                                                                &=C^{\prime}_{mk}[F^{(k)}_{J\alpha}-F^{(m)}_{J\alpha}] 
      \end{aligned}
      \end{cases}

where
:math:`C^{\prime}_{mk}=\frac{d C_{mk}(\epsilon_{mk})}{d\epsilon_{mk}}`
and :math:`F^{(k,m)}_{J\alpha}` is the :math:`\alpha` component of the
total configurational force over particle :math:`J` in the chemical
state :math:`(k,m)`. Similarly to
eq. :eq:`EevbPsimat_eq`,
eq. :eq:`Fevb2_eq` can be expressed as a matrix
multiplication

.. math::
   :label: FevbPsimat_eq

   \label{eq:FevbPsimat}
      F_{J\alpha}^{\text{EVB}}=-\sum_{m,k=1}^{N_F} \tilde{\Psi}^{(m)}_{\text{EVB}} \left(\frac{\partial H^{mk}_{\text{EVB}}}{\partial R_{J\alpha}}\right) \Psi^{(k)}_{\text{EVB}}.


The above equations define the standard EVB force field (EVB-FF). Even though the EVB formalism was first developed to compute molecular systems, EVB is also applicable to extended systems, customarily modelled using the supercell approximation and periodic boundary conditions (PBCs). However, the application of the EVB method to NPT ensembles requires the computation of the EVB stress tensor, which cannot be derived using the standard formulation :cite:`scivetti-evb`. To circumvent this limitation, we propose to make use of the well-known relation between the configurational energy and the configurational stress tensor :cite:`essmann-95a`

.. math::
   :label: stress-def1_eq

   \label{eq:stress-def1}
       \frac{\partial E^{(k)}_{c}}{\partial h_{\alpha\beta}}=-V\sum_{\gamma=x,y,z}\sigma_{\alpha\gamma}^{c(k)}h^{-1}_{\beta\gamma}

where :math:`h` is the set of lattice vectors of the supercell with
volume :math:`V`\ =det(\ :math:`h`). Multiplying to the left by
:math:`h_{\nu\beta}` and summing over :math:`\beta` we obtain the
inverse relation to eq. :eq:`stress-def1_eq`

.. math::
   :label: stress_def2_eq

   \label{eq:stress-def2}
       \sigma_{\alpha\beta}^{c(k)}=-\frac{1}{V}\sum_{\gamma=x,y,z}h_{\beta\gamma}\frac{\partial E^{(k)}_{c}}{\partial h_{\alpha\gamma}}

which can be used to define the EVB stress tensor

.. math::
   :label: stress-def3_eq

   \label{eq:stress-def3}
       \sigma_{\alpha\beta}^{\text{EVB}}=-\frac{1}{V}\sum_{\gamma=x,y,z}h_{\beta\gamma}\frac{\partial E_{\text{EVB}}}{\partial h_{\alpha\gamma}}.

Similar to the definition of the EVB force, we evaluate
:math:`\partial E_{\text{EVB}}/\partial h_{\alpha\gamma}` using the
eq. :eq:`EevbPsi_eq` and the Hellman-Feynman theorem

.. math::
   :label: stress-EVB_eq

   \label{eq:stress-EVB}
       \frac{\partial E_{\text{EVB}}}{\partial h_{\alpha\beta}}=\big\langle \Psi_{\text{EVB}}\big| \frac{\partial \hat{H}_{\text{EVB}}}{\partial h_{\alpha\beta}}\big| \Psi_{\text{EVB}}\big \rangle.

The matrix components of the operator
:math:`\frac{\partial \hat{H}_{\text{EVB}}}{\partial_{h_{\alpha\beta}}}`
follow from the definition of the EVB matrix :eq:`evbmatrix_eq` and the use of relation :eq:`stress-def1_eq`

.. math::
   :label: stress-EVB-mat_eq

   \label{eq:stress-EVB-mat}
       \frac{\partial H^{mk}_{\text{EVB}}}{\partial h_{\alpha\beta}}=\begin{cases}
         \frac{\partial E_{c}^{(m)}}{\partial h_{\alpha\beta}}=-V\sum_{\gamma}\sigma_{\alpha\gamma}^{c(m)}h^{-1}_{\beta\gamma} \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, m=k  \\
       \\
       \begin{aligned}
       \frac{d C_{mk}}{\partial h_{\alpha\beta}}&= \frac{d C_{mk}(\epsilon_{mk})}{d \epsilon_{mk}}\frac{\partial \epsilon_{mk}}{\partial h_{\alpha\beta}}
         \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \,\,\,\,\,\,\,\, \,\,\,\,\,\,\,\, \,\,\,\,\,\,\,\,    m \ne k
       \\ 
                                                                     &= \frac{d C_{mk}(\epsilon_{mk})}{d\epsilon_{mk}}\left[\frac{\partial E_{c}^{(m)}}{\partial h_{\alpha\beta}}-\frac{\partial E_{c}^{(k)}}{\partial h_{\alpha\beta}} \right]\
       \\
                                                                     &=-VC^{\prime}_{mk}\sum_{\gamma}[\sigma_{\alpha\gamma}^{c(m)}-\sigma_{\alpha\gamma}^{c(k)}] h^{-1}_{\beta\gamma}.\\
       \end{aligned}
       \end{cases} \nonumber

Finally, the EVB stress tensor of
eq. :eq:`stress-def3_eq` can be expressed as a
matrix multiplication

.. math::
   :label: stress-EVB-ab_eq

   \label{eq:stress-EVB-ab}
       \sigma_{\alpha\beta}^{\text{EVB}}=-\frac{1}{V}\sum_{\gamma=x,y,z}h_{\beta\gamma}\sum_{m,k=1}^{N_F} \tilde{\Psi}^{(m)}_{\text{EVB}} \left(\frac{\partial H^{mk}_{\text{EVB}}}{\partial h_{\alpha\beta}}\right) \Psi^{(k)}_{\text{EVB}}.

These expressions provide an alternative to compute the stress tensor
:math:`\sigma^{\text{EVB}}` from the configurational stress tensors of
each non-reactive FF, :math:`\sigma_{\alpha\gamma}^{c(k)}`. It is
important to note that the presented scheme to compute
:math:`\sigma^{\text{EVB}}` can only be derived if one uses functional
forms for :math:`C_{mk}` that depend on the energy differences
:math:`\epsilon_{mk}`, for which one can evaluate
:math:`\frac{\partial E_{c}^{(m)}}{\partial h_{\alpha\beta}}-\frac{\partial E_{c}^{(m)}}{\partial h_{\alpha\beta}}`
and use relation :eq:`stress-def1_eq` with the
computed configurational stress tensor for each chemical state. In
contrast, if the choice was to use coupling terms that do not depend on
:math:`\epsilon_{mk}` but other degrees of freedom such as spatial
coordinates, it is not clear how to derive an expression for
:math:`\sigma^{\text{EVB}}`. Similarly to the stress tensor, the
inability to compute individual contributions of the EVB force
:cite:`scivetti-evb` prevents the evaluation of the virial
using the standard formulation, and the usual decomposition of the
virial depending of the type of interaction under consideration. Within
the presented formalism, we compute the virial
:math:`\mathcal{V}_{\text{EVB}}` from
:math:`\sigma_{\alpha\beta}^{\text{EVB}}` as follows

.. math::
   :label: virial-total_eq

   \label{eq:virial-total}
       \mathcal{V}_{\text{EVB}}=-\sum_{\alpha=x,y,z} \sigma_{\alpha\alpha}^{\text{EVB}}.

The instantaneous total stress tensor, :math:`\sigma^{T}`, is given by
the following general expression

.. math::
   :label: stress-total_eq

   \label{eq:stress-total}
       \sigma^{T}=\sigma^{\text{kin}}+\sigma^{\text{EVB}}+\sigma^{\text{RB}}+\sigma^{\text{bc}}

where :math:`\sigma^{\text{kin}}`, :math:`\sigma^{\text{RB}}` and
:math:`\sigma^{\text{bc}}` are the contributions to the stress tensor
from the kinetic energy, rigid bodies (RB) and bond constraints (bc),
respectively. The EVB method only accounts for the configurational
interactions, as described. The kinetic stress tensor is computed as
usual from the instantaneous velocities of the particles. For a particle
that is part of a rigid body, the only possible interactions are
intermolecular non-bonded interactions (such as coulombic and van der
Waals interactions) with other neighboring particles that are not part
of the same rigid body. Following the computation of the EVB forces via
eq. :eq:`Fevb2_eq`, the contribution to the stress from
the rigid bodies is analogously to eq.
:eq:`rb-stress_eq`

.. math::
   :label: stress-RG_eq

   \label{eq:stress-RG}
       \sigma_{\alpha\beta}^{\text{RB}}=\sum_{\mathcal{B}=1}^{N_{\text{RB}}}\sum_{I=1}^{\eta_{\mathcal{B}}} {F}_{I_{\mathcal{B}},\alpha}^{\text{EVB}} d_{I_{\mathcal{B}},\beta}

where :math:`\vec{F}_{I_{\mathcal{B}}}` is the total force over particle
:math:`I` of rigid body :math:`\mathcal{B}` and
:math:`\vec{d}_{I_{\mathcal{B}}}` the vector distance from atom
:math:`I_{\mathcal{B}}` to the center of mass of the rigid body
:math:`\mathcal{B}`. In the above expression, index :math:`\mathcal{B}`
runs over all the rigid bodies. Each rigid body is composed of
:math:`\eta_{\mathcal{B}}` particles. Since, by definition, the topology
of rigid bodies remain unaltered during the simulation, the use of RBs
within in the present framework is meaningful only to model the
environment interacting reactive EVB site. A common example is the use
of rigidly constrained water molecules to model a solution.
Contributions to the stress tensor from bond constraints,
:math:`\sigma_{\alpha\beta}^{\text{bc}}`, are obtained using the
SHAKE/RATTLE algorithm (sec. :ref:`shake-rattle`) during
the course of the simulation. This algorithm is independent of the EVB
formalism, and corrects for the dynamics of the constrained particles.
Finally, frozen particles do not contributed to the stress tensor and
are not considered in the formalism. It is important to note that the
topology defined via the setting of RBs, frozen atoms and bond
constraints must be the consistent for all the coupled FFs, as they
impose well defined conditions for the dynamics. For example, if a group
of atoms form a rigid body, they must remain a rigid body independently
of chemical state under consideration.

.. _evb-calibrate:

Calibrating EVB force fields
----------------------------

The quality of EVB for the description of reactive processes depends on
the choice for the coupling terms :math:`C_{mk}`, particularly to
reproduce accurate interactions at the intermediate region between
chemical states :math:`m` and :math:`k` where the change of chemistry
occurs. For the implementation of the EVB method in DL_POLY_4, we have
used functional forms :math:`C_{mk}` that depend on the energy
differences :math:`\epsilon_{mk}=E^{(m)}_{c}-E^{(k)}_{c}` to compute the
stress tensor as described in Sec. :ref:`evbMethod`. We have
implemented two functional forms for the coupling terms. One is just
setting the coupling term to be a constant:

.. math::
   :label: coupl-const_eq

   \label{eq:coupl-const}
   C_{mk}(\epsilon_{mk})=\mathcal{A}_{1,mk}

and the other possibility is to use Gaussian type of function,

.. math::
   :label: coupl-gauss_eq

   \label{eq:coupl-gauss}
   C_{mk}(\epsilon_{mk})=\mathcal{A}_{1,mk} \, \, e^{-\left( \frac{\epsilon_{mk}-\mathcal{A}_{2,mk}} {\mathcal{A}_{3,mk}}  \right)^2 }+\mathcal{A}_{4,mk}.

To determine the parameters for the coupling terms, it is necessary to
consider a path that connects the reference geometries for states
:math:`m` and :math:`k`. A convenient path is the minimum energy path
(MPE) at zero-temperature, :math:`\zeta_{mk}`, obtained either via
Density Functional Theory (DFT) or quantum chemistry (QC) methods to
reproduce the change of chemistry between the states. The corresponding
energy profile for this trajectory, :math:`\tilde{E}_{\zeta_{mk}}`, is
used as a reference, and the aim is to fit the coupling parameters such
that :math:`E_{EVB}` coincides with :math:`\tilde{E}_{\zeta_{mk}}` along
:math:`\zeta_{mk}`. If we consider another state :math:`l`, for example,
it is expected that along :math:`\zeta_{mk}` the values for
:math:`E^{(l)}_{c}` will be exceedingly large in comparison with
:math:`E^{(m)}_{c}` and :math:`E^{(k)}_{c}`
(:math:`|\epsilon_{lk}|\gg 1` and :math:`|\epsilon_{lm}|\gg 1`), from
which :math:`C_{ml}(\epsilon_{ml})\approx  \mathcal{A}_{4,ml}` and
:math:`C_{kl}(\epsilon_{kl}) \approx \mathcal{A}_{4,kl}`. One can
initially set :math:`\mathcal{A}_{4,kl}=\mathcal{A}_{4,ml}=0` for all
:math:`l\ne m, k` and the coupling term :math:`C_{ml}` is computed as
follows

.. math::
   :label: coupl-neb_eq

   \label{eq:coupl-neb}
   C^{2}_{mk}(\epsilon_{mk})=\left[ \tilde{E}_{\zeta_{mk}}-E^{(m)}_{c,\zeta_{mk}} \right] \left[ \tilde{E}_{\zeta_{mk}}-E^{(k)}_{c,\zeta_{mk}} \right]

where :math:`E^{(m)}_{c,\zeta_{mk}}` and :math:`E^{(k)}_{c,\zeta_{mk}}`
are the conformational energies for states :math:`m` and :math:`k` along
:math:`\zeta_{mk}`, while :math:`\epsilon_{mk}` is in turn a implicit
function of :math:`\zeta_{mk}`

.. math::
   :label:coupl-EG_eq

   \label{eq:coupl-EG}
   \epsilon_{mk}(\zeta_{mk})=E^{(m)}_{c,\zeta_{mk}}-E^{(k)}_{c,\zeta_{mk}}.

To find the parameters for the coupling :math:`C_{mk}`, one has to plot
the values obtained from eq. :eq:`coupl-neb_eq` as a
function of :math:`\epsilon_{mk}` and fit the parameters using the
functional form of eqs. :eq:`coupl-gauss_eq` or
:eq:`coupl-const_eq`, depending on the user’s
choice. This procedure is enough when coupling two FFs. For more than
two fields, however, we have assumed :math:`\mathcal{A}_{4,ln}=0` for
:math:`l\ne n \ne m,k`. Thus, in order to fit the parameters for the
rest of the coupling terms of the EVB matrix, one should consider all
the possible remaining MEPs between the states. For the pair
:math:`l,p`, for example, one can proceed in a similar way by setting
all :math:`\mathcal{A}_{4}` elements to zero, but this time
:math:`C_{mk}` will not be necessarily zero. Depending on the number of
coupled FFs, different but more complicated expressions like eq.
:eq:`coupl-neb_eq` can be derived. Details are beyond
the scope of this chapter. The procedure to fit the coupling terms
necessarily requires the use of force-fields i) consistent with the
level of theory that is used to compute the explicit electronic problem
for the reaction and ii) accurate enough far from the reference geometry
for which they were fitted. Ultimately, meeting these requirements is a
non-trivial challenge, particularly for large systems. We refer the user
to Ref. :cite:`scivetti-evb` (and references therein) for a
more detail discussion of the available strategies to calibrate EVB
potentials. Finally, depending on the non-reactive FF and the result
from a DFT/QC simulation, one may want to shift the configuration
energies :math:`E^{(m)}_{c}` by :math:`\Delta E^{(m)}_{shift}`. This is
particularly convenient to correct the relative energy between the
involved chemical states. We have implemented this feature as input
parameters in the SETEVB file (see Sec.
:ref:`evb-users`.


.. _implement:

Computational implementation
----------------------------

In the standard format, DL_POLY_4 reads the initial coordinates,
velocities and forces from the CONFIG file. Each particle is labelled
according to its specification in the FIELD file, which contains the
information of the FF type and parameters for the interactions between
the particles. Settings for the MD simulation are specified in the
CONTROL file. Initially, the code was modified to allow i) reading
multiple (:math:`N_F`) CONFIG and FIELD files, ii) allocating arrays
of dimension :math:`N_F` for the relevant quantities, iii) checking
consistency of specification between all force fields and initial
coordinates (including any possible constraint such as rigid bodies),
iv) reading EVB settings such as coupling terms and v) preventing the
execution if there are MD or FF options that are not consistent with a
EVB simulation. With regards to this last point, not all type of
interactions in the energy decomposition of
eq. :eq:`decomp-ene` are suitable to describe
reactive interactions. For example, three-body, four-body, Tersoff and
metallic interactions are, by construction, not designed to account
for different chemical states. Thus, such interactions should only be
used to model the surrounding atomic environment interacting with the
EVB site. Regarding the EVB method in itself, modifications to the
code required to allow for the computations of energies, forces,
stress tensor and virials for each of the :math:`N_F` force-fields
separately. From the computed configurational energy of each FF and
the choice of the functional forms for the coupling terms, the EVB
matrix (eq. :eq:`evbmatrix_eq`) is built and diagonalized,
and the lowest eigenvalue and the corresponding vector are assigned to
:math:`E_{EVB}` and :math:`\Psi_{EVB}`, respectively. Matrix
:eq:`gradevb_eq` is computed for each particle’s
Cartesian components and the resulting EVB force is obtained via the
matrix multiplication of eq. :eq:`FevbPsimat_eq`.
From the stress tensors computed for each FF, matrix :eq:`stress-EVB-mat_eq`
is built for all the
:math:`\alpha\beta` terms and the :math:`\alpha\beta` component of the
EVB stress tensor obtained via
eq. :eq:`stress-EVB-ab_eq`, and the total
virial from eq. :eq:`virial-total_eq`. Such EVB
calculations are conducted for each time step taking advantage of the
domain decomposition as implemented in DL_POLY_4.
All the :math:`N_F` force fields are computed in a loop architecture,
i.e. one after the other, before being coupled via the EVB method.
This means that all the available processors are used to compute each
force-field, in contrast to the alternative strategy of dividing
processors for each force field. For extended systems, this choice is
convenient given the relative high computational cost of the long
range Coulombic part in comparison with all the other contributions to
the configurational energy. This loop structure increases the
computational time by a multiplicative factor of approximately
:math:`N_F` with respect to the required time to compute only a single
force field.

.. _evb-users:

Setting EVB calculations
------------------------

Setting input files and parameters for the EVB simulation of :math:`N_F`
coupled FFs in DL_POLY_4 requires of:

-  the CONTROL file with the directive :math:`evb\,\,\,\, N_F`

-  :math:`N_F` CONFIG files with the same ionic coordinates. The
   labelling of each atom in the CONFIG file must be consistent its
   FIELD file.

-  :math:`N_F` FIELD files with the interaction parameters to describe
   each of the coupled chemical states. Very important: the FFs
   descriptors for the reactive part of the potential must be specified
   before the descriptors for the non-reactive part.

-  the SETEVB file.

To avoid problems, users are advised to check consistency between CONFIG
and FIELD files for each of the chemical states separately, as for any
standard simulation with DL_POLY_4. It is important to remark that the
numbering and coordinates for all the atoms should be same of all CONFIG
files, and all CONFIG files must have the same number of atoms. For
example, atom 1 with tag A in CONFIG (labelling consistent with FIELD
file) should be also atom 1 in CONFIG2, even though it might have a
different tag B (labelling assigned in FIELD2). The file SETEVB is
compulsory for EVB simulations. For a EVB site described by :math:`N_F`
fields, the SETEVB file must contain all the settings specified via the
structure details in Table :numref:`(%s)<setevb_table>`. The
definition of the :math:`N_F` values of :math:`evbtypemols` in the
SETEVB file requires of particular care. As described in table
:numref:`(%s)<setevb_table>`, these values indicate how many of
the first defined type-of-molecules for each FIELD files are used to
describe the EVB reactive site. To further clarify on this statement,
let us consider a single EVB reactive unit interacting with non-reactive
water molecules. Such a reactive unit is described by two-coupled FFs.
In the chemical state 1, the reactive site is a single fragment
described by the first type-of-molecule in the FIELD file, while the
second type-of-molecule describes each of the surrounding water
molecules. In the chemical state 2, the reactive site is composed of two
molecular fragments, described by the first two type-of-molecules in the
FIELD2 file, while now the third type-of-molecule describes the
surrounding water. Consequently, Molecular types is set to 2 and 3 for
files FIELD and FIELD2, respectively, and the specification in the
SETEVB must be: :math:`evbtypemols\,\,\,\,  1\,\,\,\,  2`. The
definition of constraints is only valid for atoms that are not part of
the reactive EVB site. In addition, constraints must be kept consistent
between FIELD files. For example, if a bond-constraint is set in the
FIELD file for atoms :math:`X` and :math:`Y`, this bond-constraint
should also be defined for the other FIELD files. Similarly with
ridig-bodies, tethers, core-shells and frozen atoms. This requirement is
crucial to ensure correctness in the dynamics of the system as forces
over constrained atoms must be corrected to comply with the constraint.
In case there is an inconsistency found, the code will abort the
execution. For the non-reactive part of the system (non-EVB atoms), it
is also important to make sure that the specification for labels, mass
and charges for all non-EVB atoms is the same for all FIELD files.
Likewise, all intermolecular (Tersoff, metallic, three-body, four-body),
intramolecular (bond, angle, dihedral and inversion) and vdW
interactions between these non-EVB atoms must be the same for all FIELD
files. If any of these requirements is not fulfilled, DL_POLY_4 aborts
the execution and print an error message that (hopefully) will guide the
user to identify and fix the inconsistency. Finally, the EVB
implementation offers the possibility to restart the simulation, as
:math:`N_F` REVCON files are written. Analogous to the standard restart
calculation, the user must copy the REVIVE file to REVOLD, while each
REVCON (REVCON2, ...., REVCON\ :math:`N_F`) file must be copied to the
corresponding CONFIG (CONFIG2, ...., CONFIG\ :math:`N_F`) file. To
restart, the user must add the word :math:`restart` in CONTROL file.
Additional points for further consideration:

-  all FIELD files must have the same units

-  Replay calculations are not allowed for EVB

-  Simulations with four-body interactions are prevented

-  external electric and magnetic fields are not possible within the EVB
   formalism (see ref. :cite:`scivetti-evb`)

.. _setevb_table:

.. list-table::
   :header-rows: 1

   *  -  Setting 
      -  Desciption 
   *  -  :math:`evbtypemols`
      -  (Compulsory) Indicates how many of the firsttype-of-molecules specified in each of the :math:`N_F` FIELD files areused to describe the EVB reactive site. See :ref:`evb-users`
   *  -  :math:`evbcoupl` 
      -  (Compulsory) Specifies the information for couplingparameters. Since :math:`C_{mk}=C_{km}` only :math:`N_F(N_F-1)/2` of these lines are needed. If the specification for any pair is repeatedthe simulation is stopped. The syntax for the specification is asfollows, depending if one sets a coupling term to be a constant (:math:`const`) or use the Gaussian functional form (:math:`gauss`): 
            | :math:`evbcoul~~~m~~~k~~~const~~~\mathcal{A}_{1,mk}`
            | :math:`evbcoul~~~m~~~k~~~gauss~~~\mathcal{A}_{1,mk}~~\mathcal{A}_{2,mk}~~\mathcal{A}_{3,mk}~~\mathcal{A}_{4,mk}`
         The oder for :math:`m` and :math:`k` is irrelevant. Execution will stop if: 
            -  :math:`evbcoul` is misspelled 
            -  :math:`m=k`, :math:`m`, :math:`k<1` or :math:`m`, :math:`k>N_{F}` 
            -  Input type different from :math:`const` or :math:`gauss`.
            -  missing :math:`\mathcal{A}_{mk}` parameters 
            -  the specification of any pair is repeated 
   *  -  :math:`evbshift` 
      -  (Compulsory) specifies the energy shift for a given FF. The syntax for this is specification as follows:
            | :math:`evbshift~~~m~~~\Delta E^{(m)}_{shift}` (in units of the FIELD file)  
         Execution will stop if:
            -  :math:`evbshift` is misspelled
            -  :math:`m<1` or :math:`m>N_{F}`
            -  missing :math:`\Delta E^{(m)}_{shift}` parameters.
            -  the specification for a given FF is repeated.
   *  -  :math:`evbpop`
      - (Optional) If present, the :math:`N_F` computed values of :math:`|\Psi^{(k)}_{\text{EVB}}\big|^{2}` are printed (index :math:`(k)` for all FFs) at each time step in file POPEVB (only after equilibration). POPEVB is not overwritten upon :math:`restart`.


As an illustrative example of the SETEVB file, we consider the case of a
single reactive malonaldehyde molecule in non-reactive water.

=========== = === =================================            
evbtypemols 1 1  
evbcoupl    1 2   const 49.0 # In units of kcal/mol
evbshift    1 0.0 # In units of kcal/mol
evbshift    2 0.0 # In units of kcal/mol
=========== = === =================================

In this case, we have two possible conformations for the molecule, each
conformation described by a different FF (see section 6 of ref.
:cite:`scivetti-evb`). In contrast, the FF for the water
molecules is the same independently of the malonaldehyde conformation.
As per guidance above, the FIELD files must first specify the FF
descriptors for the malonaldehyde molecule, followed by the FF
descriptor for the surrounding water. For both FIELDS we have 2 types of
molecules, the malonaldehyde molecule for the first type and :math:`N`
water molecules for the second type. Therefore, directive
:math:`evbcoupl` must be set to :math:`1\,\,\,\,\, 1`. For the EVB
coupling, we must specify the :math:`evbcoul` with the involved fields
(1 and 2 in this case), the type of coupling (:math:`const`) and the
parameter (49.0) in units of the FIELD files. Finally, in the present
case, both conformations for malonaldehyde are energetically equivalent
and the energy shift :math:`evb` must be the same for both directives
:math:`evbshift`. Depending on the used potentials and the system to be
computed, one may want to introduce an asymmetry in the FFs.
