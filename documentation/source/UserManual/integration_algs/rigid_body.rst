.. _rigid:

Rigid Bodies and Rotational Integration Algorithms
==================================================

Description of Rigid Body Units
-------------------------------

A rigid body unit is a collection of point entities whose local geometry
is time invariant. One way to enforce this in a simulation is to impose
a sufficient number of bond constraints between the atoms in the unit.
However, in many cases this may be either problematic or impossible.
Examples in which it is impossible to specify sufficient bond
:index:`constraints<constraints;bond>` are

#. linear molecules with more than 2 atoms (e.g. CO\ :math:`_2`)

#. planar molecules with more than three atoms (e.g. benzene).

Even when the structure *can* be defined by bond constraints the network
of bonds produced may be problematic. Normally, they make the iterative
SHAKE/RATTLE procedure slow, particularly if a ring of 
:index:`constraints<constraints;bond>` is
involved (as occurs when one defines water as a constrained triangle).
It is also possible, inadvertently, to over constrain a molecule (e.g.
by defining a methane tetrahedron to have 10 rather than 9 bond
constraints) in which case the SHAKE/RATTLE procedure will become
unstable. In addition, massless sites (e.g. charge sites) cannot be
included in a simple constraint approach making modelling with
potentials such as TIP4P water impossible.

All these problems may be circumvented by defining :index:`rigid body` units, the
dynamics of which may be described in terms of the translational motion
of the centre of mass (COM) and rotation about the COM. To do this we
need to define the appropriate variables describing the position,
orientation and inertia of a :index:`rigid body`, and the rigid body equations of
motion [1]_.

The mass of a rigid unit :math:`M` is the sum of the atomic masses in
that unit:

.. math:: M = \sum_{j=1}^{N_{sites}} m_{j}~~, \label{netmass}
   :label: netmass_eq

where :math:`m_{j}` is the mass of an atom and the sum includes all
sites (:math:`N_{sites}`) in the body. The position of the rigid unit is
defined as the location of its centre of mass :math:`\underline{R}`:

.. math:: \underline{R} = {1 \over M}\sum_{j=1}^{N_{sites}} m_{j} \underline{r}_{j}~~,

where :math:`\underline{r}_{j}` is the position vector of atom :math:`j`. The
rigid body translational velocity :math:`\underline{V}` is defined by:

.. math:: \underline{V} = {1 \over M}\sum_{j=1}^{N_{sites}} m_{j} \underline{v}_{j}

and its :index:`angular momentum` :math:`\underline{J}` can then be defined by the
expression:

.. math:: \underline{J} = \sum_{j=1}^{N_{sites}} m_{j} \left( \underline{d}_{j} \times \left[ \underline{v}_{j} - \underline{V} \right] \right)~~,

where :math:`\underline{v}_{j}` is the velocity of atom :math:`j` and
:math:`\underline{d}_{j}` is the displacement vector of the atom :math:`j`
from the COM, is given by:

.. math:: \underline{d}_{j}=\underline{r}_{j}-\underline{R}~~.

The net translational force :math:`\underline{F}` acting on the :index:`rigid body`
unit is the vector sum of the forces acting on the atoms of the body:

.. math:: \underline{F} = \sum_{j=1}^{N_{sites}} \underline{f}_{j} \label{netforce}
   :label: netforce_eq

and the torque vector :math:`\underline{\tau}` acting on the body in the
universal frame of reference is given by:

.. math:: \underline{\tau} = \sum_{j=1}^{N_{sites}} \underline{d}_{j} \times \underline{f}_{j}~~,\label{nettorque}
   :label: nettorque_eq

where :math:`\underline{f}_{j}` is the force on a rigid unit site.

A rigid body also has associated with it a rotational inertia matrix
:math:`\underline{\underline{\mathbf{I}}}`, whose components are given by:

.. math::

   I^{\alpha\beta}=\sum_{j=1}^{N_{sites}} m_{j} (d_{j}^{2}
   \delta_{\alpha \beta}-d_{j}^{\alpha} r_{j}^{\beta})~~,

and COM stress and virial respectively written down as:

.. math::
   :label: rb-stress_eq

   \begin{aligned}
   \label{eq:rb-stress}
   \sigma^{\alpha\beta}=&\sum_{j=1}^{N_{sites}} d_{j}^{\alpha} f_{j}^{\beta} \nonumber \\
   {\cal W} =& -\sum_{j=1}^{N_{sites}} \underline{d}_{j} \cdot \underline{f}_{j}~~,\end{aligned}

where :math:`\underline{d}_{j}` is the displacement vector of the atom
:math:`j` from the COM, and is given by:

.. math:: \underline{d}_{j}=\underline{r}_{j}-\underline{R}~~.

The rigid body angular velocity :math:`\underline{\omega}` is the right dot
product of the inverse of the moment inertia, :math:`\underline{\underline{\mathbf{I}}}`, and the
angular momentum, :math:`\underline{J}`,:

.. math:: \underline{\omega} =  \underline{\underline{\mathbf{I}}}^{-1} \cdot \underline{J}~~.

It is common practice in the treatment of rigid body motion to define
the position :math:`\underline{R}` of the body in a universal frame of
reference (the so called laboratory or inertial frame), but to describe
the moment of inertia tensor in a frame of reference that is localised
in the rigid body and changes as the rigid body rotates. Thus the local
body frame is taken to be that in which the rotational inertia tensor
:math:`\hat{\underline{\underline{\mathbf{I}}}}` is diagonal and the components satisfy
:math:`I_{xx} \ge I_{yy} \ge I_{zz}`. In this local frame (the so called
*Principal Frame*) the inertia tensor is therefore constant.

The orientation of the local body frame with respect to the space fixed
frame is described via a four dimensional unit vector, the quaternion:

.. math:: \underline{q} = [q_0,q_1,q_2,q_3]^T~~,

and the rotational matrix :math:`\underline{\underline{\mathbf{R}}}` to transform from the local
body frame to the space fixed frame is the unitary matrix:

.. math::

   \underline{\underline{\mathbf{R}}} =
   \begin{pmatrix}
    q_0^2+q_1^2-q_2^2-q_3^2 & 2~(q_1~q_2-q_0~q_3)     & 2~(q_1~q_3+q_0~q_2)     \cr
    2~(q_1~q_2+q_0~q_3)     & q_0^2-q_1^2+q_2^2-q_3^2 & 2~(q_2~q_3-q_0~q_1)     \cr
    2~(q_1~q_3-q_0~q_2)     & 2~(q_2~q_3+q_0~q_1)     & q_0^2-q_1^2-q_2^2+q_3^2 \cr
   \end{pmatrix}

so that if :math:`\hat{\underline{d}}_{j}` is the position of an atom in the
local body frame (with respect to its COM), its position in the
universal frame (w.r.t. its COM) is given by:

.. math:: \underline{d}_{j} = \underline{\underline{\mathbf{R}}} \cdot \hat{\underline{d}}_{j}~~.

With these variables defined we can now consider the equations of motion
for the rigid body unit.

Integration of the Rigid Body Equations of Motion
-------------------------------------------------

The equations of translational motion of a rigid body are the same as
those describing the motion of a single atom, except that the force is
the total force acting on the rigid body i.e. :math:`\underline{F}` in
equation :eq:`netforce_eq` and the mass is the total mass of
the rigid body unit i.e. :math:`M` in
equation :eq:`netmass_eq`. These equations can be integrated
by the standard Verlet VV algorithm described in the previous sections.
Thus we need only consider the rotational motion here.

.. index::
   single: equations of motion;rigid body 

The rotational equation of motion for a rigid body relates the torque to
the change in angular momentum:

.. math:: \underline{\tau} = {d \over dt}\underline{J} ={d \over dt}\left(\underline{\underline{\mathbf{I}}} \cdot \underline{\omega}\right)~~.\label{torquemotion}
   :label: torquemotion

In a thermostat it can be written as:

.. math:: \dot{\underline{J}_{i}} = \underline{\tau}_{i} - \underline{\omega}_{i} \times \underline{J}_{i} + \frac{\chi}{q_{mass}} \underline{J}_{i}~~,

where :math:`i` is the index of the rigid body, :math:`\chi` and
:math:`{q_{mass}}` are the thermostat friction coefficient and
mass [2]_. In the local frame of the rigid body and without the
thermostat term, these simplify to the Euler’s 
:index:`equations<equations of motion;Euler>`

.. math::
   :label: euler_eq

   \begin{aligned}
   \dot{\hat{\omega}}_{x} =&
   {\hat{\tau}_{x} \over \hat{I}_{xx}}+(\hat{I}_{yy}-\hat{I}_{zz})~\hat{\omega}_{y}~\hat{\omega}_{z} \nonumber \\
   \dot{\hat{\omega}}_{y} =&
   {\hat{\tau}_{y}\over \hat{I}_{yy}}+(\hat{I}_{zz}-\hat{I}_{xx})~\hat{\omega}_{z}~\hat{\omega}_{z} \label{euler} \\
   \dot{\hat{\omega}}_{z} =&
   {\hat{\tau}_{z}\over \hat{I}_{zz}}+(\hat{I}_{xx}-\hat{I}_{yy})~\hat{\omega}_{x}~\hat{\omega}_{y}~~. \nonumber\end{aligned}

The vectors :math:`\hat{\underline{\tau}}` and :math:`\hat{\underline{\omega}}` are
the torque and angular velocity acting on the body transformed to the
local body frame. Integration of :math:`\hat{\underline{\omega}}` is
complicated by the fact that as the rigid body rotates, so does the
local reference frame. So it is necessary to integrate
equations :eq:`euler_eq` simultaneously with an integration of
the quaternions describing the orientation of the rigid body. The
equation describing this is:

.. math::

   \begin{pmatrix}
   \dot{q}_0 \cr \dot{q}_1 \cr \dot{q}_2 \cr \dot{q}_3 \cr
   \end{pmatrix}= {1 \over 2}
   \begin{pmatrix}
   q_0 & -q_1 & -q_2 & -q_3 \cr
   q_1 & ~q_0 & -q_3 & ~q_2 \cr
   q_2 & ~q_3 & ~q_0 & -q_1 \cr
   q_3 & -q_2 & ~q_1 & ~q_0
   \end{pmatrix}
   \begin{pmatrix}
   0 \cr \hat{\omega}_{x} \cr \hat{\omega}_{y} \cr \hat{\omega}_{z}
   \end{pmatrix}~~.


.. index:: 
   single: algorithm;FIQA 
   single: algorithm;NOSQUISH

Rotational motion in DL_POLY_4is handled by two different methods. For
the LFV implementation, the Fincham Implicit Quaternion Algorithm (FIQA)
is used :cite:`fincham-92a`. The VV implementation uses the
NOSQUISH algorithm of Miller *et al.* :cite:`miller-02a`.
The implementation NOSQUSH is coded in ``no_squish`` both contained
within ``quaternion_container``.

The LFV implementation begins by integrating the angular velocity
equation in the local frame:

.. math::

   \hat{\underline{\omega}}(t+{ \Delta t\over 2}) = \hat{\underline{\omega}}(t-{ \Delta t\over 2}) +
    \Delta t \; \hat{\underline{\underline{\mathbf{I}}}}^{-1} \cdot \dot{\hat{\underline{\omega}}}(t)~~.

The new :index:`quaternions` are found using the :index:`FIQA<algorithm;FIQA>` algorithm. In this
algorithm the new quaternions are found by solving the implicit
equation:

.. math::

   \underline{q}(t+\Delta t) = \underline{q}(t) + {\Delta t\over 2}
   \left( \underline{\underline{\mathbf{Q}}}~[\underline{q}(t)] \cdot \hat{\underline{w}}(t) +
   \underline{\underline{\mathbf{Q}}}~[\underline{q}(t+\Delta t)] \cdot \hat{\underline{w}}(t+\Delta t)\right)~~,

where :math:`\hat{\underline{w}} = [ 0,\hat{\underline{\omega}}]^T` and
:math:`\underline{\underline{\mathbf{Q}}}[\underline{q}]` is:

.. math::

   \underline{\underline{\mathbf{Q}}} = {1 \over 2}
   \begin{pmatrix}
   q_0 & -q_1 & -q_2 & -q_3 \cr
   q_1 & ~q_0 & -q_3 & ~q_2 \cr
   q_2 & ~q_3 & ~q_0 & -q_1 \cr
   q_3 & -q_2 & ~q_1 & ~q_0
   \end{pmatrix}~~.

The above equation is solved iteratively with

.. math:: {\underline{q}}(t+\Delta t) = {\underline{q}}(t) + \Delta t~\underline{\underline{\mathbf{Q}}}[{\underline{q}}(t)] \cdot \hat{\underline{w}}(t)

as the first guess. Typically, no more than 3 or 4 iterations are needed
for convergence. At each step the normalisation constraint:

.. math:: \| \underline{q}(t+\Delta t) \| = 1

is imposed.

While all the above is enough to build LFV implementations, the VV
implementations, based on the :index:`NOSQUISH<algorithm;NOSQUISH>` algorithm of Miller *et al.*
:cite:`miller-02a`, also require treatment of the quaternion
momenta as defined by:

.. math::

   \begin{pmatrix}
   p_0 \cr p_1 \cr p_2 \cr p_3 \cr
   \end{pmatrix}
   = 2
   \begin{pmatrix}
   q_0 & -q_1 & -q_2 & -q_3 \cr
   q_1 & ~q_0 & -q_3 & ~q_2 \cr
   q_2 & ~q_3 & ~q_0 & -q_1 \cr
   q_3 & -q_2 & ~q_1 & ~q_0
   \end{pmatrix}
   \begin{pmatrix}
   0 \cr \hat{I}_{xx}~\hat{\omega}_{x} \cr \hat{I}_{yy}~\hat{\omega}_{y} \cr \hat{I}_{zz}~\hat{\omega}_{z}
   \end{pmatrix}~~,

and quaternion torques as defined by:

.. math::

   \begin{pmatrix}
   \Upsilon_0 \cr \Upsilon_1 \cr \Upsilon_2 \cr \Upsilon_3 \cr
   \end{pmatrix}
   = 2
   \begin{pmatrix}
   q_0 & -q_1 & -q_2 & -q_3 \cr
   q_1 & ~q_0 & -q_3 & ~q_2 \cr
   q_2 & ~q_3 & ~q_0 & -q_1 \cr
   q_3 & -q_2 & q_1 & q_0
   \end{pmatrix}
   \begin{pmatrix}
   0 \cr \hat{\tau}_{x} \cr \hat{\tau}_{y} \cr \hat{\tau}_{z}
   \end{pmatrix}~~.

It should be noted that vectors :math:`\underline{p}` and
:math:`\underline{\Upsilon}` are 4-component vectors. The quaternion momenta
are first updated a half-step using the formula:

.. math:: \underline{p}(t+{\Delta t\over 2}) \leftarrow \underline{p}(t)+{\Delta t \over 2} \underline{\Upsilon}(t)~~.

Next a sequence of operations is applied to the quaternions and the
quaternion momenta in the order:

.. math::
   :label: ns1_eq

   e^{i{\cal L}_{3}(\delta t/2)}~e^{i{\cal L}_{2}(\delta t/2)}~e^{i{\cal L}_{1}(\delta t)}~e^{i{\cal L}_{2}(\delta t/2)}~e^{i{\cal L}_{3}(\delta t/2)}~~,
   \label{ns1}

which preserves the symplecticness of the operations (see reference
:cite:`martyna-96a`). Note that :math:`\delta t` is some
submultiple of :math:`\Delta t`. (In DL_POLY_4 the default is
:math:`\Delta t=10
\delta t`.) The operators themselves are of the following kind:

.. math::

   \begin{aligned}
   e^{i{\cal L} (\delta t)}~\underline{q}=&\cos(\zeta_{k} \delta t)~\underline{q}+\sin(\zeta_{k} \delta t)~P_{k}~\underline{q} \nonumber \\
   e^{i{\cal L} (\delta t)}~\underline{p}=&\cos(\zeta_{k} \delta t)~\underline{p}+\sin(\zeta_{k} \delta t)~P_{k}~\underline{p}~~,\end{aligned}

where :math:`P_{k}` is a permutation operator with :math:`k=0,\ldots,3`
with the following properties:

.. math::
   :label: ns2_eq

   \begin{aligned}
   P_0~\underline{q}=&\{~q_0,~q_1,~q_2,~q_3\} \nonumber \\
   P_1~\underline{q}=&\{-q_1,~q_0,~q_3,-q_2\} \label{ns2} \\
   P_2~\underline{q}=&\{-q_2,-q_3,~q_0,~q_1\} \nonumber \\
   P_3~\underline{q}=&\{-q_3,~q_2,-q_1,~q_0\}~~, \nonumber\end{aligned}

and the angular velocity :math:`\zeta_{k}` is defined as:

.. math:: \zeta_{k}={1 \over 4 I_{k}}\underline{p}^{T} P_{k}~\underline{q}~~.

Equations :eq:`ns1_eq`-\ :eq:`ns2_eq`) represent the heart of
the NOSQUISH algorithm and are repeatedly applied (10 times in ). The
final result is the quaternion updated to the full timestep value i.e.
:math:`\underline{q}(t+\Delta t)`. These equations form part of the first
stage of the VV algorithm (VV1).

In the second stage of the VV algorithm (VV2), new torques are used to
update the quaternion momenta to a full timestep:

.. math:: \underline{p}(t+\Delta t) \leftarrow \underline{p}(t+{\Delta t \over 2})+{\Delta t \over 2} \underline{\Upsilon}(t+\Delta t)~~.

Thermostats and Barostats coupling to the Rigid Body Equations of Motion
------------------------------------------------------------------------

In the presence of rigid bodies in the atomic system the system’s
instantaneous pressure, equation :eq:`prs_inst_eq`:

.. math::

   {\cal P}(t) = \frac{\left[ 2 E_{kin}(t) -
   {\cal W}_{\rm atomic}(t) - {\cal W}_{\rm COM}(t) -
   {\cal W}_{\rm constrain}(t - \Delta t) -
   {\cal W}_{\rm PMF}(t - \Delta t) \right]} {3 V(t)}

and stress, equation :eq:`str_inst_eq`:

.. math::

   \underline{\underline{\mathbf{\sigma}}}(t) = \underline{\underline{\mathbf{\sigma}}}_{kin}(t) +
   \underline{\underline{\mathbf{\sigma}}}_{\rm atomic}(t) + \underline{\underline{\mathbf{\sigma}}}_{\rm COM}(t) +
   \underline{\underline{\mathbf{\sigma}}}_{\rm constrain}(t - \Delta t) + \underline{\underline{\mathbf{\sigma}}}_{\rm PMF}(t - \Delta t)

are augmented to include the RBs’ COM virial and stress contributions.

.. note::
   
   The kinetic energy and stress in the above also include
   the contributions of the RBs’ COMs kinetic energy and stress!


.. index::
   single: barostat 
   single: thermostat

In DL_POLY_4 all degrees of freedom, translational and rotational, are
considered equal and thus treated in the same manner in all available
thermostats. Similarly, in the same spirit of equi-partitioning, all
translational degrees of freedom, the free particles’ ones and the RBs’
COMs ones are considered equal and thus treated in the same manner in
all available barostats. Based on these considerations, it is
straightforward to couple the rigid body equations of motion to a
thermostat and/or barostat. The thermostat is coupled to both the
translational and rotational degrees of freedom and thus the
translational and rotational velocities (momenta) are thermostated in
the same operational manner as the purely atomic ones. The barostat,
however, is coupled only to the translational degrees of freedom and
does not contribute to the rotational motion of the RBs, thus only the
RBs’ COMs positions and momenta are subjected to the same barostat
driven algorithmic operations as those of the free particles’ positions
and momenta. Therefore, if we notion the change of the system’s degrees
of freedom as:

.. math:: f \rightarrow F = f + f^{RB \mathbf{ (tra)}} + f^{RB \mathbf{ (rot)}}

then all equations of motion defining the ensembles as described in this
chapter are subject to the following notional changes in order to
include the RB contributions:

.. math::

   \begin{aligned}
   \sigma(f) \rightarrow& ~~~~~\sigma(F) ~~~~~~~~~~~~~~~~~~~~= ~~~~~~\sigma(f + f^{RB \mathbf{ (tra)}} + f^{RB \mathbf{ (rot)}}) \nonumber \\
   {\cal H}(f) \rightarrow& ~~~~{\cal H}(F) ~~~~~~~~~~~~~~~~~~~~= ~~~~~{\cal H}(f + f^{RB \mathbf{ (tra)}} + f^{RB \mathbf{ (rot)}}) \nonumber \\
   p_{mass}(f) \rightarrow& p_{mass}(F - f^{RB \mathbf{ (rot)}}) = ~~~~~p_{mass}(f + f^{RB \mathbf{ (tra)}}) \\
   \eta(f) \rightarrow& ~~~~~\eta(F - f^{RB \mathbf{ (rot)}})~~ = ~~~~~~\eta(f + f^{RB \mathbf{ (tra)}}) \nonumber \\
   \underline{\underline{\mathbf{\eta}}}(f) \rightarrow& ~~~~~\underline{\underline{\mathbf{\eta}}}(F - f^{RB \mathbf{ (rot)}})~~ = ~~~~~~\underline{\underline{\mathbf{\eta}}}(f + f^{RB \mathbf{ (tra)}})~~, \nonumber\end{aligned}

where :math:`f` refers to the degrees of freedom in the system (see
equation :eq:`freedom_eq`), :math:`\sigma` is the system
target energy (see equation :eq:`sigma_eq`), :math:`{\cal H}`
is the conserved quantity of the ensemble (if there is such defined),
:math:`E_{kin}` (!includes RB COM kinetic energy too) and
:math:`E_{rot}` are respectively the kinetic and rotational energies of
the system, :math:`p_{mass}` is the barostat mass, and :math:`\eta` and
:math:`\underline{\underline{\mathbf{\eta}}}` are the barostat friction coefficient or matrix of
coefficients respectively.

There are two slight technicalities with the Evans and Andersen
ensembles that are worth mentioning.

Since both the translational and rotational velocities contribute
towards temperature, equation :eq:`Evans_eq`, showing the
derivation of the thermostat friction in the Evans ensemble by imposing
a Gaussian constraint on the system’s instantaneous temperature, changes
to:

.. math::

   \begin{aligned}
   \frac{d}{dt} {\cal T} = 0 ~~~~~~\propto~~~~~~ \frac{d}{dt} \left( \frac{1}{2} \sum_{i}^{FP} m_{i} \underline{v}_{i}^{2} + \frac{1}{2} \sum_{j}^{RB} M_{j} \underline{V}_{j}^{2} +
   \frac{1}{2} \sum_{j}^{RB} \hat{\underline{\omega}}_{j}^{T} \cdot \hat{\underline{\underline{\mathbf{I}}}}_{j} \cdot \hat{\underline{\omega}}_{j} \right) = 0 \nonumber \\
   \left\{ \sum_{i}^{FP} \underline{v}_{i}(t) \cdot \underline{f}_{i}(t) + \sum_{j}^{RB} \underline{V}_{j}(t) \cdot \underline{F}_{j}(t) +
   \sum_{j}^{RB} \hat{\underline{\omega}}_{j}(t) \cdot \hat{\underline{\tau}}(t) \right\} - ~~~~~~~~~~~~~~~~~~~~~~~~~  \nonumber \\
   ~~~~~~~~~~~~~~~~~~~\chi (t) \left\{ \sum_{i}^{FP} m_{i} \underline{v}_{i}^{2}(t) +
   \sum_{j}^{RB} M_{j} \underline{V}_{j}^{2}(t) + \sum_{j}^{RB} \hat{\underline{\omega}}_{j}^{T}(t) \cdot \hat{\underline{\underline{\mathbf{I}}}}_{j} \cdot \hat{\underline{\omega}}_{j}(t) \right\} = 0 \\
   \chi (t) = \frac {\left\{ \sum_{i}^{FP} \underline{v}_{i}(t) \cdot \underline{f}_{i}(t) + \sum_{j}^{RB} \underline{V}_{j}(t) \cdot \underline{F}_{j}(t) +
   \sum_{j}^{RB} \hat{\underline{\omega}}_{j}(t) \cdot \hat{\underline{\tau}}(t) \right\}} {2\left[ E_{kin}(t) + E_{rot}(t) \right]}~~, \nonumber\end{aligned}

where where :math:`\cal T` is the instantaneous temperature defined in
equation :eq:`tinst_eq` and :math:`E_{kin}` in the final
expression contains both the kinetic contribution form the free
particles and the RBs’ COMs.

In the case of the Andersen ensemble, if a Poisson selected particle
constitutes a RB then the whole RB is Poisson selected. Poisson selected
RBs’ translational and angular velocities together with Poisson selected
FPs’ velocities sample the same Gaussian distribution isokinetically
(Boltzmann distribution), where the isokineticity to target temperature
is dependent upon the total of the Poisson selected FPs’ and RBs’
degrees of freedom.


.. [1]
   An alternative approach is to define “basic” and “secondary”
   particles. The basic particles are the minimum number needed to
   define a local body axis system. The remaining particle positions are
   expressed in terms of the COM and the basic particles. Ordinary bond
   :index:`constraints<constraints;bond>` can then be applied to the 
   basic particles provided the
   forces and torques arising from the secondary particles are
   transferred to the basic particles in a physically meaningful way.

.. [2]
   It is worth noting that in DL_POLY_4all degrees of freedom,
   translational (both the free particles’ ones and the RBs’ COMs ones)
   and rotational, are considered equal and thus treated in the same
   manner in all available thermostats!