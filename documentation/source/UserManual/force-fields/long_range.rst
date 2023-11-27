.. _coulomb:

.. index:: single: potential;electrostatics

Long Ranged Electrostatic (coulombic) Potentials
================================================

.. index:: single: potential;electrostatics

DL_POLY_4 incorporates several techniques for dealing with long-ranged
electrostatic potentials [1]_. These are as follows:

#. Direct Coulomb sum

#. Force-shifted Coulomb sum

#. Coulomb sum with distance dependent dielectric

#. Reaction field

#. Smoothed Particle Mesh Ewald (SPME)

.. index:: single: polarisation;shell model

All of these can be used in conjunction with the shell model technique
used to account for ions polarisation.

The SPME technique is restricted to periodic systems only. (Users must
exercise care when using pseudo-periodic :index:`boundary conditions`.) The other
techniques can be used with either periodic or non-periodic systems
safely, although in the case of the :index:`direct Coulomb sum` there are likely
to be problems with convergence.

DL_POLY_4 will correctly handle the electrostatics of both molecular and
atomic species. However, it is assumed that the system is electrically
neutral. A warning message is printed if the system is found to be
charged, but otherwise the simulation proceeds as normal.

.. note::
   
   DL_POLY_4 does not use the basic Ewald method, which is an
   option in , on account of it being too slow for large scale systems. The
   SPME method is the standard Ewald method in .

Default (Point Charges) Electrostatics
--------------------------------------

.. index:: single: direct Coulomb sum

Direct Coulomb Sum
~~~~~~~~~~~~~~~~~~

Use of the direct Coulomb sum is sometimes necessary for accurate
simulation of isolated (non-periodic) systems. It is *not* recommended
for periodic systems.

The interaction potential for two charged ions is

.. math:: U(r_{ij}) = \frac{1}{4\pi\epsilon_{0}\epsilon}\frac{q_{i}q_{j}}{r_{ij}}~~,

with :math:`q_{\ell}` the charge on an atom labelled :math:`\ell`, and
:math:`r_{ij}` the magnitude of the separation vector
:math:`\underline{r}_{ij}=\underline{r}_{j}-\underline{r}_{i}` .

The force on an atom :math:`j` derived from this force is

.. math:: \underline{f}_{j} = \frac{1}{4\pi\epsilon_{0}\epsilon} \frac{q_{i}q_{j}}{r_{ij}^{3}} \underline{r}_{ij}~~,

with the force on atom :math:`i` the negative of this.

The contribution to the atomic virial is

.. math:: {\cal W} = -\frac{1}{4\pi\epsilon_{0}\epsilon} \frac{q_{i}q_{j}}{r_{ij}}~~,

which is simply the negative of the potential term.

The contribution to be added to the atomic :index:`stress tensor` is

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha} f_{j}^{\beta}~~,

where :math:`\alpha,\beta` are :math:`x,y,z` components. The atomic
:index:`stress tensor` is symmetric.

In DL_POLY_4 these forces are handled by the subroutine ``coul_cp_forces``.

.. index:: force-shifted Coulomb sum

Force-Shifted Coulomb Sum
~~~~~~~~~~~~~~~~~~~~~~~~~

This form of the Coulomb sum has the advantage that it drastically
reduces the range of electrostatic interactions, without giving rise to
a violent step in the potential energy at the cutoff. Its main use is
for preliminary preparation of systems and it is not recommended for
realistic models.

The form of the simple truncated and shifted potential function is

.. math::

   U(r_{ij}) = \frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon}
   \left\{\frac{1}{r_{ij}} - \frac{1}{r_{\rm cut}}\right\}~~,

with :math:`q_{\ell}` the charge on an atom labelled :math:`\ell`,
:math:`r_{\rm cut}` the cutoff radius and :math:`r_{ij}` the magnitude
of the separation vector :math:`\underline{r}_{ij}=\underline{r}_{j}-\underline{r}_{i}` .

A further refinement of this approach is to truncate the :math:`1/r`
potential at :math:`r_{\rm cut}` and add a linear term to the potential
in order to make both the energy and the force zero at the cutoff. This
removes the heating effects that arise from the discontinuity in the
forces at the cutoff in the simple truncated and shifted potential (the
formula above). (The physics of this potential, however, is little
better. It is only recommended for very crude structure optimizations.)

The force-shifted potential is thus

.. math::

   U(r_{ij}) = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon}
   \left[ \left\{\frac{1}{r_{ij}}+\frac{1}{r_{\rm cut}^{2}}~r_{ij}\right\} -
   \left\{\frac{1}{r_{\rm cut}}+\frac{1}{r_{\rm cut}^{2}}~r_{\rm cut}\right\} \right]
   = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon}
   \left[ \frac{1}{r_{ij}} + \frac{r_{ij}}{r_{\rm cut}^{2}} - \frac{2}{r_{\rm cut}} \right]~~,

with the force on an atom :math:`j` given by

.. math::

   \underline{f}_{j} = \frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon}
   \left[ \frac{1} {r_{ij}^{3}} - \frac{1}{r_{ij}r_{\rm cut}^{2}} \right] \underline{r}_{ij}~~,

with the force on atom :math:`i` the negative of this.

The force-shifted Coulomb potential can be elegantly extended to emulate
long-range ordering by including distance depending damping function
:math:`{\rm erfc}(\alpha~r_{ij})` (identical to that seen in the
real-space portion of the Ewald sum) and thus mirror the effective
charge screening :cite:`fennell-06a` as shown below

.. math::

   \begin{aligned}
   U(r_{ij}) = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon} \left[
   \left\{ \frac{{\rm erfc}(\alpha~r_{ij})}{r_{ij}} + \left(\frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}}+
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}}\right)~r_{ij}\right\} - \phantom{xx} \right. \nonumber \\
   \left. \left\{ \frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}} + \left(\frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}}+
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}}\right)~r_{\rm cut}\right\} \right]~~,\end{aligned}

with the force on an atom :math:`j` given by

.. math::

   \begin{aligned}
   \underline{f}_{j} = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon} \left[ \left( \frac{{\rm erfc}(\alpha~r_{ij})}{r_{ij}^{2}} +
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{ij}^{2})}{r_{ij}} \right) - \phantom{xxxxx} \right. \nonumber \\
   \left. \left( \frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}} +
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}} \right) \right]~\frac{\underline{r}_{ij}}{r_{ij}}~~,
   \end{aligned}

with the force on atom :math:`i` the negative of this.

It is worth noting that, as discussed in :cite:`fennell-06a`
and references therein, this is only an approximation of the Ewald sum
and its accuracy and effectiveness become better when the cutoff is
large (\ :math:`>` 10 preferably 12 Å).

The contribution to the atomic virial is

.. math:: {\cal W} = -\underline{r}_{ij} \cdot \underline{f}_{j}~~,

which is *not* the negative of the potential term in this case.

The contribution to be added to the atomic :index:`stress tensor` is given by

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha} f_{j}^{\beta}~~,

where :math:`\alpha,\beta` are :math:`x,y,z` components. The atomic
stress tensor is symmetric.

In DL_POLY_4 these forces are handled by the routine ``coul_fscp_forces``.

.. index:: single: distance dependant dielectric

Coulomb Sum with Distance Dependent Dielectric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This potential attempts to address the difficulties of applying the
:index:`direct Coulomb sum`, without the brutal truncation of the previous case.
It hinges on the assumption that the electrostatic forces are
effectively ‘screened’ in real systems - an effect which is approximated
by introducing a dielectric term that increases with distance.

The interatomic potential for two charged ions is

.. math:: U(r_{ij}) = \frac{1}{4\pi\epsilon_{0}\epsilon(r_{ij})} \frac{q_{i}q_{j}}{r_{ij}}~~,

with :math:`q_{\ell}` the charge on an atom labelled :math:`\ell`, and
:math:`r_{ij}` the magnitude of the separation vector
:math:`\underline{r}_{ij}=\underline{r}_{j}-\underline{r}_{i}` . :math:`\epsilon(r)` is
the distance :index:`dependent dielectric` function. In DL_POLY_4 it is assumed
that this function has the form

.. math:: \epsilon(r)~=~\epsilon~r~~,

where :math:`\epsilon` is a constant. Inclusion of this term effectively
accelerates the rate of convergence of the Coulomb sum.

The force on an atom :math:`j` derived from this potential is

.. math:: \underline{f}_{j} = \frac{1}{2\pi\epsilon_{0}\epsilon} \frac{q_{i}q_{j}}{r_{ij}^{4}} \underline{r}_{ij}~~,

with the force on atom :math:`i` the negative of this.

The contribution to the atomic virial is

.. math:: {\cal W} = -\underline{r}_{ij} \cdot \underline{f}_{j}~~,

which is :math:`-2` times the potential term.

The contribution to be added to the atomic :index:`stress tensor` is given by

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha} f_{j}^{\beta}~~,

where :math:`\alpha,\beta` are :math:`x,y,z` components. The atomic
:index:`stress tensor` is symmetric.

In DL_POLY_4 these forces are handled by the routine ``coul_dddp_forces``.

.. index:: single: reaction field

Reaction Field
~~~~~~~~~~~~~~

In the :index:`reaction field` method it is assumed that any given molecule is
surrounded by a spherical cavity of finite radius within which the
electrostatic interactions are calculated explicitly. Outside the cavity
the system is treated as a :index:`dielectric<distance dependant dielectric>` continuum. The occurrence of any
net dipole within the cavity induces a polarisation in the dielectric,
which in turn interacts with the given molecule. The model allows the
replacement of the infinite Coulomb sum by a finite sum plus the
reaction field.

The reaction field model coded into DL_POLY_4 is the implementation of
Neumann based on charge-charge interactions
:cite:`neumann-85a`. In this model, the total coulombic
potential is given by

.. math::

   U_{c} = \frac{1}{4\pi\epsilon_{0}\epsilon} \sum_{j<n} q_{j}q_{n} \left[
   \frac{1}{r_{nj}} + \frac{B_{0}r_{nj}^{2}}{2 R_{c}^{3}} \right]~~,

where the second term on the right is the reaction field correction to
the explicit sum, with :math:`R_{c}` the radius of the cavity. The
constant :math:`B_{0}` is defined as

.. math:: B_{0} = \frac{2(\epsilon_{1}-1)}{(2\epsilon_{1}+1)}~~,

with :math:`\epsilon_{1}` the dielectric constant outside the cavity.
The effective pair potential is therefore

.. math::

   U(r_{ij}) = \frac{1}{4\pi\epsilon_{0}\epsilon} q_{i}q_{j} \left[
   \frac{1}{r_{ij}} + \frac{B_{0}r_{ij}^{2}}{2 R_{c}^{3}} \right]~~.

This expression unfortunately leads to large fluctuations in the system
coulombic energy, due to the large ‘step’ in the function at the cavity
boundary. In DL_POLY_4 this is countered by subtracting the value of the
potential at the cavity boundary from each pair contribution. The term
subtracted is

.. math::

   \frac{1}{4\pi\epsilon_{0}\epsilon} \frac{q_{i}q_{j}}{R_{c}} \left[
   1+\frac{B_{0}}{2} \right]~~.

The effective pair force on an atom :math:`j` arising from another atom
:math:`n` within the cavity is given by

.. math::

   \underline{f}_{j}=\frac{q_{i}q_{j}}{4\pi\epsilon_{0}\epsilon}\left[
   \frac{1}{r_{ij}^{3}}-\frac{B_{0}}{R_{c}^{3}}\right] \underline{r}_{ij}~~.

In DL_POLY_4 the reaction field is optionally extended to emulate
long-range ordering in a force-shifted manner by countering the reaction
term and using a distance depending damping function
:math:`{\rm erfc}(\alpha~r_{ij})` (identical to that seen in the
real-space portion of the Ewald sum) and thus mirror the effective
charge screening :cite:`fennell-06a`:

.. math::

   \begin{aligned}
   U(r_{ij}) = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon}& \left[
   \left\{ \frac{{\rm erfc}(\alpha~r_{ij})}{r_{ij}} + \left(\frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}}+
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}}\right)~r_{ij}\right\} \right. \\
   &- \left\{ \frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}} + \left(\frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}}+
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}}\right)~r_{\rm cut}\right\} \nonumber \\
   &+ \left. \frac{B_{0}(r_{ij}^{2}-r_{\rm cut}^{2})}{2 r_{\rm cut}^{3}} \right]~~, \nonumber
   \end{aligned}

with the force on an atom :math:`j` given by

.. math::

   \begin{aligned}
   \underline{f}_{j} = {q_{i} q_{j} \over 4\pi\epsilon_{0}\epsilon} \left[ \left( \frac{{\rm erfc}(\alpha~r_{ij})}{r_{ij}^{2}} +
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{ij}^{2})}{r_{ij}} \right) - \phantom{xxxxxxxxxxx} \right. \\
   \left. \left( \frac{{\rm erfc}(\alpha~r_{\rm cut})}{r_{\rm cut}^{2}} +
   \frac{2\alpha}{\sqrt{\pi}}~\frac{\exp(-\alpha^{2}~r_{\rm cut}^{2})}{r_{\rm cut}} \right) -
   \frac{B_{0}r_{ij}}{r_{\rm cut}^{3}} \right]~\frac{\underline{r}_{ij}}{r_{ij}}~~, \nonumber\end{aligned}

with the force on atom :math:`i` the negative of this.

It is worth noting that, as discussed in :cite:`fennell-06a`
and references therein, this is only an approximation of the Ewald sum
and its accuracy and effectiveness become better when the cutoff is
large (:math:`>` 10 preferably 12 Å).

The contribution of each effective pair interaction to the atomic virial
is

.. math:: {\cal W} = -\underline{r}_{ij} \cdot \underline{f}_{j}

and the contribution to the atomic :index:`stress tensor` is

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha} f_{j}^{\beta}~~,

where :math:`\alpha,\beta` are :math:`x,y,z` components. The atomic
:index:`stress tensor` is symmetric.

In DL_POLY_4 the reaction field is handled by the subroutine
``coul_rfp_forces``.

.. _SPME:

Smoothed Particle Mesh Ewald
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: single: Ewald;summation

The Ewald sum :cite:`allen-89a` is the best technique for
calculating electrostatic interactions in a periodic (or
pseudo-periodic) system.

The basic model for a neutral periodic system is a system of charged
point ions mutually interacting via the Coulomb potential. The Ewald
method makes two amendments to this simple model. Firstly, each ion is
effectively neutralised (at long-ranged) by the superposition of a
spherical :index:`Gaussian<constraints;Gaussian>` cloud of opposite charge 
centred on the ion. The
combined assembly of point ions and :index:`Gaussian<constraints;Gaussian>`
charges becomes the *Real
Space* part of the :index:`Ewald<Ewald;summation>` sum, which is now short ranged and treatable by
the methods described above
(Chapter :ref:`force-field`) [2]_. The second modification
is to superimpose a second set of Gaussian charges, this time with the
same charges as the original point ions and again centred on the point
ions (so nullifying the effect of the first set of Gaussians). The
potential due to these :index:`Gaussian<constraints;Gaussian>`s 
is obtained from Poisson’s equation and
is solved as a Fourier series in *Reciprocal Space*. The complete 
:index:`Ewald<Ewald;summation>`
sum requires an additional correction, known as the self energy
correction, which arises from a :index:`Gaussian<constraints;Gaussian>` 
acting on its own site, and is
constant. Ewald’s method, therefore, replaces a potentially infinite sum
in real space by two finite sums: one in real space and one in
reciprocal space; and the self energy correction.

For molecular systems, as opposed to systems comprised simply of point
ions, additional modifications
are necessary to correct for the excluded (intra-molecular) coulombic
interactions. In the real space sum these are simply omitted. In
reciprocal space however, the effects of individual :index:`Gaussian<constraints;Gaussian>` charges
cannot easily be extracted, and the correction is made in real space.
It amounts to removing terms corresponding to the potential energy of
an ion :math:`\ell` due to the :index:`Gaussian<constraints;Gaussian>` charge on a neighbouring
charge :math:`m` (or *vice versa*). This correction appears in the
term noting a summation over :math:`molecules` in the full :index:`Ewald<Ewald;summation>`
formula below.

The same considerations and modifications ewald_frzn_forces are taken
into account for frozen atoms, which mutual coulombic interaction must
also be excluded. This correction appears in the term noting a summation
over :math:`F^{*}` (all frozen-frozen pairs in the MD cell) in the full
:index:`Ewald<Ewald;summation>` formula below.

Note the distinction between the *error function* **erf** and the more
usual *complementary error function* **erfc** found in the real space
sums below.

The total electrostatic energy is given by the following formula:

.. math::

   \begin{aligned}
   U_{c}=&\frac{1}{2V_{o} \epsilon_{0}\epsilon} \sum_{\underline{k} \neq
   \underline{0}}^{\underline{\infty}} \frac{\exp(-k^{2}/4\alpha^{2})}{k^{2}}
   \left|\sum_{j}^{N} q_{j}\exp(-i\underline{k}\cdot\underline{r}_{j})\right|^{2} \nonumber \\
   &-\frac{1}{4\pi\epsilon_{0}\epsilon} \frac{\alpha}{\sqrt{\pi}}
   \sum_{j}^{N} q_{j}^2  \nonumber \\
   &+ \frac{1}{4\pi\epsilon_{0}\epsilon} \sum_{n<j}^{N^{*}} \frac{q_{j}q_{n}}
   {r_{nj}} {\rm erfc}(\alpha r_{nj}) \nonumber \\
   &- \frac{1}{4\pi\epsilon_{0}\epsilon}
   \sum_{molecules} \sum_{\ell\le m}^{M^{*}} q_{\ell}q_{m} \left\{ \delta_{\ell m}
   \frac{\alpha}{\sqrt{\pi}} + \frac{{\rm erf}(\alpha r_{\ell m})}{r_{\ell
   m}^{1-\delta_{\ell m}}} \right\} \\
   &- \frac{1}{4\pi\epsilon_{0}\epsilon} \sum_{\ell\le m}^{F^{*}}
   q_{\ell}q_{m} \left\{ \delta_{\ell m} \frac{\alpha}{\sqrt{\pi}} +
   \frac{{\rm erf}(\alpha r_{\ell m})}{r_{\ell m}^{1-\delta_{\ell m}}}
   \right\} \nonumber \\
   &- \frac{1}{4\pi\epsilon_{0}\epsilon} \frac{\pi}{2V_{o} \alpha^{2}}
   \left\{ \sum_{j}^{N} q_{j} \right\}^{2}~~, \nonumber
   \end{aligned}

where :math:`N` is the number of ions in the system and :math:`N^{*}`
the same number discounting any excluded (intramolecular and frozen)
interactions. :math:`M^{*}` represents the number of excluded atoms in a
given molecule. :math:`F^{*}` represents the number of frozen atoms in
the MD cell. :math:`V_{o}` is the simulation cell volume and
:math:`\underline{k}` is a reciprocal lattice vector defined by

.. math:: \underline{k} = \ell \underline{u} + m \underline{v} + n \underline{w} \label{k-vector}~~,

where :math:`\ell,m,n` are integers and :math:`\underline{u},\underline{v},\underline{w}`
are the *reciprocal space* basis vectors. Both :math:`V_{o}` and
:math:`\underline{u},\underline{v},\underline{w}` are derived from the vectors
(:math:`\underline{a},\underline{b},\underline{c}`) defining the simulation cell. Thus

.. math:: V_{o} = |\underline{a} \cdot \underline{b} \times \underline{c}|

and

.. math::

   \begin{aligned}
   \underline{u} =& 2 \pi \frac{\underline{b} \times \underline{c}}{\underline{a} \cdot \underline{b} \times \underline{c}} \nonumber \\
   \underline{v} =& 2 \pi \frac{\underline{c} \times \underline{a}}{\underline{a} \cdot \underline{b} \times\underline{c}} \\
   \underline{w} =& 2 \pi \frac{\underline{a} \times \underline{b}}{\underline{a} \cdot \underline{b} \times \underline{c}}~~. \nonumber\end{aligned}

With these definitions, the Ewald formula above is applicable to general
periodic systems. The last term in the Ewald formula above is the Fuchs
correction :cite:`fuchs-35a` for electrically non-neutral MD
cells which prevents the build-up of a charged background and the
introduction of extra pressure due to it.

In practice the convergence of the Ewald sum is controlled by three
variables: the real space cutoff :math:`r_{\rm cut}`; the convergence
parameter :math:`\alpha` and the largest reciprocal space vector
:math:`\underline{k}_{max}` used in the reciprocal space sum. These are
discussed more fully in
Section  :ref:`ewald-precision`. DL_POLY_4 can
provide estimates if requested (see CONTROL file description
:ref:`control-file`).

.. index:: single: Ewald;SPME

As its name implies the Smoothed Particle Mesh Ewald (SPME) method is a
modification of the standard Ewald method. DL_POLY_4 implements the SPME
method of Essmann *et al.* :cite:`essmann-95a`. Formally,
this method is capable of treating van der Waals forces also, but in
DL_POLY_4 it is confined to electrostatic forces only. The main
difference from the standard Ewald method is in its treatment of the
reciprocal space terms. By means of an interpolation procedure involving
(complex) B-splines, the sum in reciprocal space is represented on a
three dimensional rectangular grid. In this form the Fast Fourier
Transform (FFT) may be used to perform the primary mathematical
operation, which is a 3D convolution. The efficiency of these procedures
greatly reduces the cost of the reciprocal space sum when the range of
:math:`\underline{k}` vectors is large. The method (briefly) is as follows
(for full details see :cite:`essmann-95a`):

#. Interpolation of the :math:`\exp(-i~\underline{k}\cdot\underline{r}_{j})` terms
   (given here for one dimension):

   .. math::

      \exp(2\pi i~u_{j} k/L) \approx b(k) \sum_{\ell=-\infty}^{\infty}
      M_{n}(u_{j}-\ell)~\exp(2\pi i~k\ell/K)~~,

   in which :math:`k` is the integer index of the :math:`\underline{k}` vector
   in a principal direction, :math:`K` is the total number of grid
   points in the same direction and :math:`u_{j}` is the fractional
   coordinate of ion :math:`j` scaled by a factor :math:`K` (i.e.
   :math:`u_{j}=K s_{j}^{x}`) . **Note** that the definition of the
   B-splines implies a dependence on the integer :math:`K`, which limits
   the formally infinite sum over :math:`\ell`. The coefficients
   :math:`M_{n}(u)` are B-splines of order :math:`n` and the factor
   :math:`b(k)` is a constant computable from the formula:

   .. math::

      b(k) = \exp(2\pi i~(n-1)k/K) \left[ \sum_{\ell=0}^{n-2}
      M_{n}(\ell+1)~\exp(2\pi i~k\ell/K) \right]^{-1}~.

#. Approximation of the structure factor :math:`S(\underline{k})`:

   .. math:: S(\underline{k}) \approx b_{1}(k_{1})~b_{2}(k_{2})~b_{3}(k_{3})~Q^{\dagger}(k_{1},k_{2},k_{3})~~,

   where :math:`Q^{\dagger}(k_{1},k_{2},k_{3})` is the discrete Fourier
   transform of the *charge array* :math:`Q(\ell_{1},\ell_{2},\ell_{3})`
   defined as

   .. math::

      \begin{aligned}
      Q(\ell_{1},\ell_{2},\ell_{3})=& \sum_{j=1}^{N}q_{j}
      \sum_{n_{1},n_{2},n_{3}} M_{n}(u_{1j}-\ell_{1}-n_{1}L_{1})~\times~
      M_{n}(u_{2j}-\ell_{2}-n_{2}L_{2}) \phantom{xxxx} \nonumber \\
      & \phantom{xxxxxxxxxx}~\times~ M_{n}(u_{3j}-\ell_{3}-n_{3}L_{3})~~,\end{aligned}

   in which the sums over :math:`n_{1,2,3}` etc are required to capture
   contributions from all relevant periodic cell images (which in
   practice means the nearest images).

#. Approximating the reciprocal space energy :math:`U_{recip}`:

   .. math::

      U_{recip} = \frac{1}{2V_{o} \epsilon_{0}\epsilon} \sum_{k_{1},k_{2},k_{3}}
      G^{\dagger}(k_{1},k_{2},k_{3})~Q(k_{1},k_{2},k_{3})~~,

   where :math:`G^{\dagger}` is the discrete Fourier transform of the
   function

   .. math::

      G(k_{1},k_{2},k_{3}) = \frac{\exp(-k^{2}/4\alpha^{2})}{k^{2}}~
      B(k_{1},k_{2},k_{3})~(Q^{\dagger}(k_{1},k_{2},k_{3}))^{*}~~,

   in which :math:`(Q^{\dagger}(k_{1},k_{2},k_{3}))^{*}` is the complex
   conjugate of :math:`Q^{\dagger}(k_{1},k_{2},k_{3})` and

   .. math:: B(k_{1},k_{2},k_{3}) = |b_{1}(k_{1})|^{2}~|b_{2}(k_{2})|^{2}~|b_{3}(k_{3})|^{2}~~.

   The function :math:`G(k_{1},k_{2},k_{3})` is thus a relatively simple
   product of the Gaussian screening term appearing in the conventional
   Ewald sum, the function :math:`B(k_{1},k_{2},k_{3})` and the discrete
   Fourier transform of :math:`Q(k_{1},k_{2},k_{3})`.

#. Calculating the atomic forces, which are given formally by:

   .. math::

      f_{j}^{\alpha} = -\frac{\partial U_{recip}}{\partial
      r_{j}^{\alpha}} = -\frac{1}{V_{o} \epsilon_{0}\epsilon}
      \sum_{k_{1},k_{2},k_{3}} G^{\dagger}(k_{1},k_{2},k_{3})
      ~\frac{\partial Q(k_{1},k_{2},k_{3})}{\partial r_{j}^{\alpha}}~~.

Fortunately, due to the recursive properties of the B-splines, these
formulae are easily evaluated.

The virial and the stress tensor are calculated in the same manner as
for the conventional Ewald sum.

The DL_POLY_4 subroutines required to calculate the SPME contributions
are:

#. spme_container containing

   #. bspgen, which calculates the B-splines

   #. bspcoe, which calculates B-spline coefficients

   #. spl_cexp, which calculates the FFT and B-spline complex
      exponentials

#. parallel_fft and gpfa_module (native DL_POLY_4 subroutines that
   respect the domain decomposition concept) which calculate the 3D
   complex fast Fourier transforms

#. ewald_spme_forces, which calculates the reciprocal space
   contributions (uncorrected)

#. ewald_real_forces, which calculates the real space contributions
   (corrected)

#. ewald_excl_forces, which calculates the reciprocal space corrections
   due to the coulombic exclusions in intramolecular interactions

#. ewald_frzn_forces, which calculates the reciprocal space corrections
   due to the exclusion interactions between frozen atoms

#. two_body_forces, in which all of the above subroutines are called
   sequentially and also the Fuchs correction
   :cite:`fuchs-35a` for electrically non-neutral MD cells
   is applied if needed.

.. _mpoles:

.. index:: multipolar electrostatics
   
Multipolar Electrostatics
-------------------------

DL_POLY_4 offers advanced potential energy calculations through
multipolar electrostatics. This is an extension to the point-charge
model where the charge density of chemical species are described by
higher order point multipoles. The generic algorithms in DL_POLY_4 are
designed to allow for arbitrary order :cite:`boateng-15a`
multipoles but for practical reasons the functionality is limited to
hexadecapoles only.

Multipoles
~~~~~~~~~~

Define the multipolar operator, :math:`\hat{L}_i` as

.. math::

   \hat{L}_i= (q_i + {\mathbf{p}}_i\cdot \nabla_i + {\mathbf{Q}}_i: \nabla_i\nabla_i +
   {\mathbf{O}}_i{\vdots} \nabla_i\nabla_i\nabla_i +
   {\mathbf{H}}_i:: \nabla_i\nabla_i\nabla_i\nabla_i + \dots)~~,

where :math:`q_i`, :math:`{\mathbf{p}}_i`, :math:`{\mathbf{Q}}_i`,
:math:`{\mathbf{O}}_i`, and :math:`{\mathbf{H}}_i` are the point
charge, dipole, quadrupole, octupole, and hexadecapole tensors,
respectively of atom *i*, :math:`\nabla_i` refers to the
three-dimensional gradient with respect to the position of atom *i*
and the “dot" products stand for tensor contraction. By defining a
unidimensional vector of independent (non-degenerate) multipole
moments, :math:`\mathcal{M}_i`, for atom *i*, the corresponding
multipolar operator to an arbitrary order :math:`p` can be written in
a more compact form as

.. math::
   :label: defLi_eq

   \hat{L}_i= \sum_{||\mathbf{s}||= 0}^{p}\mathcal{M}_{i}^{\mathbf{s}}\partial_{i}^{\mathbf{s}} =
   \sum_{s_3 = 0}^{p}\sum_{s_2 = 0}^{p-s_3}\sum_{s_1=0}^{p-s_3-s_2} \mathcal{M}_{j}^{s_1 s_2 s_3}
   {\partial}_{z_i}^{s_3}{\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}~~.\label{defLi}

Here, :math:`\mathbf{s}= (s_1,s_2,s_3)` is the triplet that runs over
all independent multipoles, :math:`||\mathbf{s}||= s_1 + s_2 + s_3`,
:math:`\mathcal{M}_{i}^{\mathbf{s}}=\mathcal{M}_{i}^{s_1 s_2 s_3}` and
:math:`\partial_{i}^{\mathbf{s}} = {\partial}_{z_i}^{s_3}{\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}`
is the multidimensional derivative with respect to the position
:math:`\langle x_i, y_i, z_i \rangle` of atom *i* with orders
:math:`s_1`, :math:`s_2` and :math:`s_3` in the :math:`x`, :math:`y`
and :math:`z` directions respectively. Individual components of
:math:`\mathcal{M}` contain the sum of all degenerate original
multipole components. As an example, the octupole
:math:`\mathcal{M}^{111}`, is a sum of all six degenerate original
octupole components formed from the permutation of the triplet
:math:`\{x,y,z\}` . If the original octupole vector with degnerate
components is labelled as :math:`O'`, then
:math:`\mathcal{M}^{111}= O_{xyz}' + O_{xzy}' + O_{yxz}' + O_{yzx}' + O_{zxy}' + O_{zyx}' = 6~O_{xyz}'` .
For pair potentials it is often convenient to redefine the multipolar
operator for atom *j* in terms of the derivatives with respect to the
position of atom *i* to arrive at

.. math::
   :label: defLj_eq

   \hat{L}_{j_i}= \sum_{||\mathbf{s}||= 0}^{p}\mathcal{M}_{j}^{\mathbf{s}}\partial_{j}^{\mathbf{s}} =
   \sum_{||\mathbf{s}||= 0}^{p}(-1)^{||\mathbf{s}||}\mathcal{M}_{j}^{\mathbf{s}}\partial_{i}^{\mathbf{s}} =
   \sum_{s_3 = 0}^{p}\sum_{s_2 = 0}^{p-s_3}\sum_{s_1=0}^{p-s_3-s_2}
   (-1)^{s_1+s_2+s_3}\mathcal{M}_{j}^{s_1 s_2 s_3}
   {\partial}_{z_i}^{s_3}{\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}~~.\label{defLj}

.. _apptopairpot:

Application to Pair Potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In DL_POLY_4 for :math:`N` point-multipoles interacting via a pair
potential function :math:`\psi`, the multipolar electrostatic potential
at position :math:`\mathbf{r_i}` is computed as

.. math::
   :label: mpolpot_eq

   \phi(\mathbf{r_i}) =  \sum_{j \ne i}^{N}\hat{L}_{j_i}\psi(\mathbf{r_{ji}}) = \sum_{j \ne i}^{N} \sum_{\mathbf{s}= \mathbf{0}}^{p}
   (-1)^{||\mathbf{s}||}\mathcal{M}_{j}^{\mathbf{s}}\partial_{i}^{\mathbf{s}}\psi(r_{ij})~~,\label{eqn:mpolpot}

the electrostatic field at :math:`\mathbf{r_i}` is

.. math::
   :label: mpolefield_eq

   \mathbf{E}(\mathbf{r_{ij}}) = -\nabla_i \phi(r_{ij})= -\sum_{j \ne i}^{N} \sum_{\mathbf{s}= \mathbf{0}}^{p}
   (-1)^{||\mathbf{s}||}\mathcal{M}_{j}^{\mathbf{s}}\left[\begin{array}{l} \partial_{i}^{\mathbf{s}+ \mathbf{e}_1} \\
   \partial_{i}^{\mathbf{s}+ \mathbf{e}_2} \\
   \partial_{i}^{\mathbf{s}+ \mathbf{e}_3}
   \end{array}\right]\psi(r_{ij})~~,\label{eqn:mpolefield}

where :math:`\mathbf{e}_1=\langle1,0, 0\rangle`,
:math:`\mathbf{e}_2=\langle0,1,0\rangle`, and
:math:`\mathbf{e}_3=\langle0,0,1\rangle` and the torque
:cite:`sagui-04a` on particle :math:`i` in the
:math:`\alpha`-direction, :math:`\tau_{i,\alpha}`, is obtained as

.. math::

   \tau_{i,\alpha} = \sum_{\mathbf{s}= \mathbf{0}}^{p} \mathcal{M}_{i,\alpha}^{\mathbf{s}} \partial_{i}^{\mathbf{s}} \phi(\mathbf{r_{ij}})
   = \sum_{\mathbf{s}= \mathbf{0}}^{p} \mathcal{M}_{i,\alpha}^{\mathbf{s}} \sum_{j \ne i}^{N} \sum_{\mathbf{k}= \mathbf{0}}^{p}
   (-1)^{||\mathbf{k}||}\mathcal{M}_{j}^{\mathbf{k}}\partial_{i}^{\mathbf{s}+\mathbf{k}}\psi(r_{ij})~~,

where :math:`\mathcal{M}_{i,\alpha}` is the infinitesimal
counter-clockwise rotation of multipole vector :math:`\mathcal{M}_i`
about the :math:`\alpha`-axis. The total electrostatic potential energy
is given by

.. math::
   :label: mpolene_eq

   U = \sum_{i < j}^{N} \hat{L}_i\hat{L}_{j_i}\psi(r_{ij}) = \sum_{i < j}^{N} \sum_{\mathbf{s}= \mathbf{0}}^{p}
   (-1)^{||\mathbf{s}||}\mathcal{M}_{j}^{\mathbf{s}}\sum_{\mathbf{k}= \mathbf{0}}^{p}\mathcal{M}_{i}^{\mathbf{k}}\partial_{i}^{\mathbf{s}+\mathbf{k}}\psi(r_{ij})~~,\label{eqn:mpolene}

where :math:`\mathbf{s}+ \mathbf{k}= (s_1+k_1,s_2+k_2,s_3+k_3)` and the
force on atom :math:`i` is

.. math::
   :label: mpolforce_eq

   \mathbf{f}_i = -\nabla_i \sum_{j \ne i}^{N} \hat{L}_i\hat{L}_{j_i}\psi(r_{ij}) =
               -\sum_{j \ne i}^{N} \sum_{\mathbf{s}= \mathbf{0}}^{p}(-1)^{||\mathbf{s}||}\mathcal{M}_{j}^{\mathbf{s}}
                \sum_{\mathbf{k}= \mathbf{0}}^{p}\mathcal{M}_{i}^{\mathbf{k}}
   \left[\begin{array}{l}
   \partial_{i}^{\mathbf{s}+\mathbf{k}+ \mathbf{e}_1} \\
   \partial_{i}^{\mathbf{s}+\mathbf{k}+ \mathbf{e}_2} \\
   \partial_{i}^{\mathbf{s}+\mathbf{k}+ \mathbf{e}_3}
   \end{array}\right] \psi(r_{ij})~~.\label{eqn:mpolforce}

To implement
equations :eq:`mpolpot_eq`-\ :eq:`mpolforce_eq`
for the variety of potentials in DL_POLY_4 a number of recurrence
relations are used to compute the multi-dimensional derivatives of the
kernels corresponding to the potentials. These kernels are

.. math::

   \begin{aligned}
   \theta(|\mathbf{x}|)=& \frac{1}{{|\mathbf{x}|}^{\nu}},\mbox{\hskip 10pt}\Omega(|\mathbf{x}|) \\ 
   =& \frac{1}{2}\textrm{exp}(-\alpha^2|\mathbf{x}|^2), \mbox{\hskip 10pt }
    \psi(|\mathbf{x}|) \\
   =&\frac{\sqrt{\pi}}{2} \frac{\textrm{erfc}(\alpha|\mathbf{x}|)}{|\mathbf{x}|},
   \mbox{\hskip 10pt and \hskip 10pt }\Gamma(|\bar{\mathbf{x}}|) \\ 
   =& \frac{\sqrt{\pi}}{2} \frac{\textrm{erf}(\alpha |\mathbf{x}|)}{|\mathbf{x}|}~~;
   \end{aligned}

with

.. math::

   \begin{aligned}
   a_{\mathbf{s}}(\nu)&=\frac{\partial^{||\mathbf{s}||}\theta(|\mathbf{x}|)}{\partial_{x_1}^{s_1}\partial_{x_2}^{s_2}\partial_{x_3}^{s_3}},\mbox{\hskip 10pt}
   b_{\mathbf{s}} \\ 
   &=\frac{\partial^{||\mathbf{s}||}\Omega(|\mathbf{x}|)}{\partial_{x_1}^{s_1}\partial_{x_2}^{s_2}\partial_{x_3}^{s_3}},\mbox{\hskip 10pt}
   c_{\mathbf{s}} \\ 
   &=\frac{\partial^{||\mathbf{s}||}\psi(|\mathbf{x}|)}{\partial_{x_1}^{s_1}\partial_{x_2}^{s_2}\partial_{x_3}^{s_3}},\mbox{\hskip 10pt and \hskip 10pt}
   d_{\mathbf{s}} \\
   &=\frac{\partial^{||\mathbf{s}||}\Gamma(|\bar{\mathbf{x}}|)}{\partial_{x_1}^{s_1}\partial_{x_2}^{s_2}\partial_{x_3}^{s_3}}~~.
   \end{aligned}

The recurrence relations used in DL_POLY_4 are

.. math::
   :label: coulrecur_eq

   a_{\mathbf{s}}(\nu) = \frac{1}{|\mathbf{x}|^2}\left\{\left(\frac{2-\nu}{||\mathbf{s}||} - 2\right)
   \sum_{i=1}^{3}s_i x_i a_{\mathbf{s}-\mathbf{e}_i} + \left(\frac{2-\nu}{||\mathbf{s}||} - 1\right)
   \sum_{i=1}^{3}s_i (s_i-1) a_{\mathbf{s}-2\mathbf{e}_i} \right\}~~,\label{coulrecur}

.. math::
   :label: exprecur_eq

   b_{\mathbf{s}} = \frac{-2\alpha^2}{||\mathbf{s}||} \sum_{i=1}^{3}
   \left[ s_i x_i b_{\mathbf{s}-\mathbf{e}_i} + s_i (s_i-1) b_{\mathbf{s}-2\mathbf{e}_i} \right]~~,\label{eqn:exprecur}

.. math::
   :label: erfcrecur_eq

   c_{\mathbf{s}} = \frac{1}{|\mathbf{x}|^2} \left\{\left(\frac{1}{||\mathbf{s}||} - 2\right)
   \sum_{i=1}^{3}s_i x_i c_{\mathbf{s}-\mathbf{e}_i} + \left(\frac{1}{||\mathbf{s}||} - 1\right)
   \sum_{i=1}^{3}s_i (s_i-1) c_{\mathbf{s}-2\mathbf{e}_i} + \frac{1}{\alpha} b_{\mathbf{s}} \right\}~~,\label{eqn:erfcrecur}

and

.. math::
   :label: erfrecur_eq

   d_{\mathbf{s}} = \frac{1}{|\mathbf{x}|^2}\left\{\left(\frac{1}{||\mathbf{s}||} - 2\right)
   \sum_{i=1}^{3}s_i x_i d_{\mathbf{s}-\mathbf{e}_i} + \left(\frac{1}{||\mathbf{s}||} - 1\right)
   \sum_{i=1}^{3}s_i (s_i-1) d_{\mathbf{s}-2\mathbf{e}_i} -\frac{1}{\alpha} b_{\mathbf{s}} \right\}~~.\label{eqn:erfrecur}

Direct Coulomb Sum
~~~~~~~~~~~~~~~~~~

For two interacting ions :math:`i` and :math:`j`, the potential energy
is given as

.. math:: U(r_{ij}) = \frac{1}{4\pi\epsilon_0\epsilon}\hat{L}_i\hat{L}_{j_i}\left[\frac{1}{r_{ij}}\right]~~,

and the relevant kernel is :math:`\psi(r_{ij}) = \frac{1}{r_{ij}}` . The
derivatives for this kernel are obtained by using
equation :eq:`coulrecur_eq` with :math:`\nu = 1` . Thus,

.. math:: \partial_i^{\mathbf{s}}\psi(r_{ij}) = a_{\mathbf{s}}(1)~~.

In DL_POLY_4 the multipolar direct Coulomb sum is handled by the routine
``coul_cp_mforces``.

Force-Shifted Coulomb Sum
~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 employs two forms of the force-shifted Coulomb sum. In the
first form, the potential energy due to two interacting ions :math:`i`
and :math:`j` is

.. math:: U(r_{ij}) = \frac{1}{4\pi\epsilon_0\epsilon} \hat{L}_i\hat{L}_{j_i}\left[\frac{1}{r_{ij}}+\frac{r_{ij}}{r_{\textrm{cut}}^2}-\frac{2}{r_{\textrm{cut}}}\right]~~,

where :math:`r_{\textrm{cut}}` is the cutoff radius. The kernel is
:math:`\psi(r_{ij}) = \frac{1}{r_{ij}}+\frac{r_{ij}}{r_{\textrm{cut}}^2}-\frac{2}{r_{\textrm{cut}}}` .
The last term, :math:`\frac{2}{r_{\textrm{cut}}}`, is a constant which
has a zero derivative, hence the derivatives of the kernel are obtained
as a sum of the derivatives of the first term and second terms. Thus,

.. math:: \partial_i^{\mathbf{s}}\psi(r_{ij}) = a_{\mathbf{s}}(1)+\frac{a_{\mathbf{s}}(-1)}{r_{\textrm{cut}}^2}~~.

The potential energy due to two point-multipoles :math:`i` and :math:`j`
interacting via the second form of the force-shifted Coulomb sum is

.. math::

   \begin{aligned}
   U(r_{ij}) =& \frac{1}{4\pi\epsilon_0\epsilon}\hat{L}_i\hat{L}_{j_i}\left[ \left\{ \frac{\textrm{erfc}(\alpha \cdot r_{ij})}{r_{ij}} +
   \left( \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} + \frac{2\alpha}{\sqrt{\pi}}
   \frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}}\right)r_{ij}\right\} \right. \nonumber \\
   & \left. \phantom{xxxxxxxxx} - \left\{ \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}} +
   \left( \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} +
   \frac{2\alpha}{\sqrt{\pi}}\frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}}\right)r_{\textrm{cut}}\right\} \right]~~.
   \end{aligned}

The kernel, :math:`\psi(r_{ij})` is the terms in the square bracket but
the only terms which contribute to the derivatives are the first and
second terms which are functions of :math:`r_{ij}` . The derivative of
the first term is obtained from
equations :eq:`erfcrecur_eq` and the derivative
for :math:`r_{ij}` in the second term is given by
:math:`d_{\mathbf{s}}(-1)` . Thus,

.. math::

   D_i^{\mathbf{s}}\psi(r_{ij}) = \frac{2}{\sqrt{\pi}}c_{\mathbf{s}} + \left(\frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} +
   \frac{2\alpha}{\sqrt{\pi}} \frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}} \right) \cdot a_{\mathbf{s}}(-1)~~.

In DL_POLY_4 the multipolar force-shifted Coulomb sum is handled by the
routine ``coul_fscp_mforces``.

Coulomb Sum with Distance Dependent Dielectric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The potential energy between two interacting ions :math:`i` and
:math:`j` is

.. math:: U(r_{ij}) = \frac{1}{4\pi\epsilon_0\epsilon}\hat{L}_i\hat{L}_{j_i}\left[\frac{1}{r_{ij}^2}\right]~~,

and the kernel is :math:`\psi(r_{ij}) = \frac{1}{r_{ij}^2}` . The
derivatives for this kernel are obtained by using
equation :eq:`coulrecur_eq` with :math:`\nu = 2` . Hence,

.. math:: \partial_i^{\mathbf{s}}\psi(r_{ij}) = a_{\mathbf{s}}(2)~~.

In DL_POLY_4 the multipolar Coulomb sum with distance dependent
dielectric is handled by the routine coul_dddp_mforces.

Reaction Field
~~~~~~~~~~~~~~

DL_POLY_4 provides two forms of a multipolar reaction field potential.
In the first form, the effective pair potential energy due to two
interacting point multipoles :math:`i` and :math:`j` is given as

.. math::

   U(r_{ij}) = \frac{1}{4\pi\epsilon_0\epsilon}\hat{L}_i\hat{L}_{j_i}\left[ \frac{1}{r_{ij}} +
   \frac{B_0r_{ij}^2}{2R_c^3} - 1 - \frac{B_0}{2}\right]~~,

where

.. math:: B_0 = \frac{2(\epsilon_1 - 1)}{(2\epsilon_1 + 1)}~~,

:math:`R_c` is the radius of the spherical cavity and :math:`\epsilon_1`
is the dielectric constant outside the cavity. Again the kernel
:math:`\psi(r_{ij})` is the terms in the square bracket and only the
first and second terms contribute to its derivatives. The derivatives of
the first and second terms are given by
equation :eq:`coulrecur_eq` with :math:`\nu = 1` and
:math:`\nu = -2` respectively. Thus,

.. math:: \partial_i^{\mathbf{s}}\psi(r_{ij}) = a_{\mathbf{s}}(1)+\frac{B_0}{2R_c^3} \cdot a_{\mathbf{s}}(-2)~~.

The second form of the reaction field method is similar to that of the
force-shifted Coulomb sum. The potential energy due to interacting ions
:math:`i` and :math:`j` is

.. math::

   \begin{aligned}
   U(r_{ij}) =& \frac{1}{4\pi\epsilon_0\epsilon}\hat{L}_i\hat{L}_{j_i}\left[ \left\{ \frac{\textrm{erfc}(\alpha \cdot r_{ij})}{r_{ij}} +
   \left( \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} + \frac{2\alpha}{\sqrt{\pi}}
   \frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}}\right)r_{ij}\right\} - \phantom{xxxxxxxxx} \right. \\ \nonumber
   & \left. \phantom{xxxxxxxxx} \left\{ \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}} +
   \left( \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} +
   \frac{2\alpha}{\sqrt{\pi}}\frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}}\right)r_{\textrm{cut}}\right\} -
   \frac{B_0 r_{\textrm{cut}}^2}{2r_{\textrm{cut}}^3}+\frac{B_0r_{ij}^2}{2r_{\textrm{cut}}^3} \right]~~.
   \end{aligned}

The kernel, :math:`\psi(r_{ij})` is the terms in the square bracket and
the only terms which contribute to the derivatives are the first, second
and last terms which are functions of :math:`r_{ij}` . The derivative of
the first term is obtained from
equation :eq:`erfcrecur_eq` and the derivative for
:math:`r_{ij}` in the second term is given by :math:`a_{\mathbf{s}}(-1)`
and the derivative for :math:`r_{ij}^2` in the last term is given by
:math:`d_{\mathbf{s}}(-2)` . Thus,

.. math::

   D_i^{\mathbf{s}}\psi(r_{ij}) = \frac{2}{\sqrt{\pi}}c_{\mathbf{s}}+\left( \frac{\textrm{erfc}(\alpha \cdot r_{\textrm{cut}})}{r_{\textrm{cut}}^2} +
   \frac{2\alpha}{\sqrt{\pi}}\frac{\textrm{exp}(-\alpha^2r_{\textrm{cut}}^2)}{r_{\textrm{cut}}}\right) \cdot a_{\mathbf{s}}(-1) +
   \frac{B_0}{2r_{\textrm{cut}}^3} \cdot a_{\mathbf{s}}(-2)~~.

In DL_POLY_4 the multipolar reaction field is handled by the routine
coul_rfp_mforces.

Smoothed Particle Mesh Ewald
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 provides two different smooth particle Mesh Ewald
implementations for multipolar electrostatics. The first implementation
is for systems with charges, dipoles and quadrupoles and does not use
recurrence relations. The second implementation, which uses recurrence
relations, is more general and allows for specification of an arbitrary
order up to hexadecapoles.

When the multipolar form of SPME is employed, the total electrostatic
energy for a system on :math:`N` point ions is given as

.. math:: U_c = U_{\textrm{dir}} + U_{\textrm{rec}} - U_{\textrm{excl}} - U_{\textrm{frzn}} - U_{\textrm{self}}~~,\label{eqn:totE}
   :label: SPME_totE_eq

where

.. math::
   :label: direwald_eq

   U_{\textrm{dir}} = \sum_{i < j}^{N^*} \sum_{\mathbf{n}}^{'}\hat{L}_i\hat{L}_{j_i}
   \frac{\textrm{erfc}(\alpha \cdot |\mathbf{r_{ij}}+ \mathbf{n}|)}{4\pi\epsilon_0\epsilon|\mathbf{r_{ij}}+ \mathbf{n}|}~~,\label{direwald}

.. math::
   :label: exclewald_eq

   U_{\textrm{excl}} = \frac{1}{4\pi\epsilon_0\epsilon}\sum_{(i,j)\in M^*}
   \hat{L}_i\hat{L}_{j_i}\frac{\textrm{erf}(\alpha \cdot r_{ij})}{r_{ij}}~~,\label{exclewald}

.. math::
   :label: frznewald_eq

   U_{\textrm{frzn}} = \frac{1}{4\pi\epsilon_0\epsilon}\sum_{(i,j)\in F^*}
   \hat{L}_i\hat{L}_{j_i}\frac{\textrm{erf}(\alpha \cdot r_{ij})}{r_{ij}}~~,\label{frznewald}

.. math::
   :label: selfewald_eq

   U_{\textrm{self}} = \frac{1}{8\pi\epsilon_0\epsilon}\lim_{| \mathbf{r_i}|\to 0}\sum_{i = 1}^{N}
   \hat{L}_i\hat{L}_i\frac{\textrm{erf}(\alpha \cdot |\mathbf{r_i}|)}{|\mathbf{r_i}|}~~,\label{selfewald}

and

.. math::
   :label: recewald_eq

   U_{\textrm{rec}} = \frac{1}{2V_o\epsilon_0\epsilon}\displaystyle\sum_{\mathbf{k} \ne 0}
   \frac{\textrm{exp}(-k^2/4\alpha^2)}{k^2}\left |S(\mathbf{k})\right|^2~~,\label{recewald}

with

.. math:: S(\mathbf{k}) = \sum_{i=1}^N \hat{L}_i\textrm{exp}(\imath \mathbf{k}\cdot \mathbf{r_i})~~.\label{eqn:sfac}
   :label: sfac_eq

In the expressions above, :math:`M^*` is the set of all excluded
interactions due to intramolecular bonds in the simulation cell,
:math:`F^*` the set of frozen-frozen interactions in the simulation
cell, :math:`N^* = N - M^* - F^*`, :math:`V_o` is the volume of the
simulation cell and :math:`S(\mathbf{k})` is the structure factor.

Real Space Sum
~~~~~~~~~~~~~~

The relevant kernel for the real space from
equation :eq:`direwald_eq` is
:math:`\displaystyle \psi(r_{ij}) = \frac{\textrm{erfc}(\alpha |\mathbf{r_{ij}}+ \mathbf{n}|)}{|\mathbf{r_{ij}}+ \mathbf{n}|}` .
DL_POLY_4 uses the recurrence giving in
equation :eq:`erfcrecur_eq` to generate the
multidimensional derivatives of the kernel. Thus, the derivatives of the
kernel are computed as

.. math:: \mathbf{D}_i^{\mathbf{s}}\psi(r_{ij}) = \frac{2}{\sqrt{\pi}}c_{\mathbf{s}}~~.

In DL_POLY_4 the routine ewald_real_mforces_d computes the real space
interactions explicitly for simulations with multipoles of order 2
without using the recurrence relation. The routine ewald_real_mforces
handles the general version of up to order 4 using recurrence relations.

Excluded Sum
~~~~~~~~~~~~

The relevant kernel for the real space from
equation :eq:`exclewald_eq` is
:math:`\displaystyle \psi(r_{ij}) = \frac{\textrm{erf}(\alpha \cdot r_{ij})}{r_{ij}}` .
DL_POLY_4 uses the recurrence giving in
equation :eq:`erfrecur_eq` to generate the
multidimensional derivatives of the kernel. Thus, the derivatives of the
kernel are computed as

.. math:: \mathbf{D}_i^{\mathbf{s}}\psi(r_{ij}) = \frac{2}{\sqrt{\pi}}d_{\mathbf{s}}~~.

In DL_POLY_4 the routine ewald_excl_mforces_d computes the reciprocal
space corrections due to the exclusions between intramolecularly related
atoms explicitly for simulations with multipoles of order 2 without
using the recurrence relation. The routine ewald_excl_mforces handles
the general version of up to order 4 using recurrence relations.

Frozen Sum
~~~~~~~~~~

The relevant kernel for the real space from
equation :eq:`frznewald_eq` is
:math:`\displaystyle \psi(r_{ij}) = \frac{\textrm{erf}(\alpha \cdot r_{ij})}{r_{ij}}` .
DL_POLY_4 uses the recurrence giving in
equation :eq:`erfrecur_eq` to generate the
multidimensional derivatives of the kernel. Thus, the derivatives of the
kernel are computed as

.. math:: \mathbf{D}_i^{\mathbf{s}}\psi(r_{ij}) = \frac{2}{\sqrt{\pi}}d_{\mathbf{s}}~~.

In DL_POLY_4 the routine ewald_frzn_mforces computes computes the
reciprocal space corrections due to the exclusions between frozen atoms
generically for simulations with multipoles up to order 4 using
recurrence relations.

Self-Interaction
~~~~~~~~~~~~~~~~

DL_POLY_4 computes :math:`U_{self}` directly for interactions involving
multipoles up to order 4 using the series representation of the kernel
:math:`\displaystyle \psi(r_{ij}) = \frac{\textrm{erf}(\alpha \cdot r_i)}{r_i}` .
The self interaction is computed in ewald_real_mforces_d for simulations
with multipoles of maximum order 2. For simulations of arbitrary order,
the self-interaction is computed in the reciprocal space.

Reciprocal Space Sum
~~~~~~~~~~~~~~~~~~~~

The key idea of SPME is in approximating the structure factor, in a
uniform grid, with :math:`K_1 \times K_2 \times K_3` dimensions, that
fills the simulation cell. Define the fractional coordinates of an ion
:math:`i` as
:math:`\langle s_{i_{1}}, s_{i_{2}}, s_{i_{3}} \rangle = \langle \mathbf{a}_1^*\cdot \mathbf{r_i}, \mathbf{a}_2^*\cdot \mathbf{r_i}, \mathbf{a}_3^*\cdot \mathbf{r_i}\rangle`,
:math:`u_{\alpha_i} = K_{\alpha} \cdot s_i^{\alpha}` and :math:`M_n`
is a B-spline of order :math:`n` then the approximation of the
structure factor is given as

.. math:: S(\mathbf{k}) \approx b_1(k_1)b_2(k_2)b_3(k_3) Q^{\mathcal{F}}(k_1,k_2,k_3)~~,\label{skapprox}
   :label: skapprox_eq

where :math:`\mathbf{k}= \langle k_1, k_2, k_3 \rangle` is a
reciprocal space vector,

.. math:: b_i(k_i) = \textrm{exp}(2\pi \imath (n-1)k_i/K_i)\left[\sum_{l=0}^{n-2}M_n(l+1)\textrm{exp}(2\pi \imath k l/K_i)\right ]^{-1},

:math:`Q` is the multipolar array defined on the uniform grid and
:math:`Q^{\mathcal{F}}` its discrete Fourier transform. At position
:math:`(l_1,l_2,l_3)` on the grid, the multipolar array is defined by

.. math::
   :label: mparray1_eq

   \begin{aligned}
   Q(l_1,l_2,l_3)=\sum_{i=1}^{N}\hat{L}_i\sum_{n_1,n_2,n_3} &M_n(u_{1_i}-l_1-n_1 K_1) \times
   M_n(u_{2_i}-l_2-n_2 K_2) \\
   &\times M_n(u_{3_i}-l_3-n_3 K_3)~~,
   \end{aligned}
   \label{marray1}

where, :math:`u_{\alpha_i}-l_{\alpha}-n_{\alpha}K_{\alpha}` are
evaluation points of the B-spline on the grid that spans the
fundamental cell and the periodic images. Then from
equation :eq:`defLi_eq` and considering only the fundamental
cell, the multipolar array can be written explicitly as

.. math::
   :label: mparray2_eq

   \begin{aligned}
   Q(l_1,l_2,l_3) = \sum_{i=1}^{N}\sum_{s_3 = 0}^{p}\sum_{s_2 = 0}^{p-s_3}
   \sum_{s_1=0}^{p-s_3-s_2}\mathcal{M}_{i}^{s_1 s_2 s_3} {\partial}_{z_i}^{s_3}
   {\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}
   \left\{ M_n(u_{1_i}-l_1) M_n(u_{2_i}-l_2) M_n(u_{3_i}-l_3)\right\}~~.
   \end{aligned}
   \label{marray2}

To compute the arbitrary order multidimensional derivatives of the
product of three b-splines in equation :eq:`mparray2_eq`,
DL_POLY_4 uses the closed form formula:

.. math::
   :label: dprodmn_eq

   \begin{aligned}
    & {\partial}_{z_i}^{s_3}{\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}
        \left\{ M_n(u_{1_i}-l_1) M_n(u_{2_i}-l_2) M_n(u_{3_i}-l_3)\right\} = \nonumber \\
    & \displaystyle \sum_{k_3 = 0}^{s_3}\left(K_1 a_{13}^{*}\right)^{k_3}
        \binom{s_3}{k_3} \sum_{k_2=0}^{s_2}\left(K_1 a_{12}^{*}\right)^{k_2}
        \binom{s_2}{k_2}\sum_{k_1=0}^{s_1}\left(K_1 a_{11}^{*}\right)^{k_1}
        \binom{s_1}{k_1} {\partial}_{u_{1_i}}^{||\mathbf{k}||}M_n(u_{1_i}-l_1) \times \\
    & \displaystyle \sum_{j_3=0}^{s_3-k_3}\left(K_2 a_{23}^{*}\right)^{j_3}
        \left(K_3 a_{33}^{*}\right)^{s_3-k_3-j_3}\binom{s_3-k_3}{j_3}
        \sum_{j_2=0}^{s_2-k_2}\left(K_2 a_{22}^{*}\right)^{j_2}
        \left(K_3 a_{32}^{*}\right)^{s_2-k_2-j_2} \binom{s_2-k_2}{j_2} \times \nonumber \\
    & \displaystyle \sum_{j_1=0}^{s_1-k_1}\left(K_2 a_{21}^{*}\right)^{j_1}
        \left(K_3 a_{31}^{*}\right)^{s_1-k_1-j_1} \binom{s_1-k_1}{j_1}
        {\partial}_{u_{2_i}}^{||\mathbf{j}||}M_n(u_{2_i}-l_2){\partial}_{u_{3_i}}^{||\mathbf{s}-\mathbf{k}-\mathbf{j}||}M_n(u_{3_i}-l_3)~~,\nonumber \label{dprodmn}
   \end{aligned}

where
:math:`\mathbf{a}_{1}^{*} = \langle a_{11}^{*},a_{12}^{*},a_{13}^{*}\rangle`,
:math:`\mathbf{a}_{2}^{*} = \langle a_{21}^{*},a_{22}^{*},a_{23}^{*}\rangle`,
and
:math:`\mathbf{a}_{3}^{*} = \langle a_{31}^{*},a_{32}^{*},a_{33}^{*}\rangle`
are the reciprocal space basis vectors and :math:`K_1`, :math:`K_2`, and
:math:`K_3`, the maximum number of grid points in the fundamental cell
in the :math:`x`, :math:`y`, and :math:`z` directions respectively. For
an orthogonal box, where

.. math:: a_{12}^{*}=a_{13}^{*}=a_{21}^{*}=a_{23}^{*}=a_{31}^{*}=a_{32}^{*}=0~~,

DL_POLY_4 uses the simplification of
equationv :eq:`dprodmn_eq` to

.. math::
   :label: dprodmnsimple_eq

   \begin{aligned}
   \label{dprodmnsimple}
    & {\partial}_{z_i}^{s_3}{\partial}_{y_i}^{s_2}{\partial}_{x_i}^{s_1}
        \left\{ M_n(u_{1_i}-l_1) M_n(u_{2_i}-l_2) M_n(u_{3_i}-l_3)\right\}= \\
    & \left(K_1 a_{11}^{*}\right)^{s_1}\left(K_2 a_{22}^{*}\right)^{s_2}
        \left(K_3 a_{33}^{*}\right)^{s_3}{\partial}_{u_{1_i}}^{s_1}M_n(u_{1_i}-l_1)
        {\partial}_{u_{2_i}}^{s_2}M_n(u_{2_i}-l_2){\partial}_{u_{3_i}}^{s_3}M_n(u_{3_i}-l_3)~~. \nonumber
   \end{aligned}

The formulas in equations :eq:`dprodmn_eq` and :eq:`dprodmnsimple_eq` require derivatives of a
b-spline. To compute an arbitrary :math:`p_\textrm{th}` order derivative
of a b-spline of order :math:`n`, :math:`M_n`, at an arbitrary grid
point :math:`j`, DL_POLY_4 uses the closed form formula

.. math::
   :label: dmnj2_eq

   \frac{d^p}{d u^{p}}M_n(u_j) = \sum_{t=\textrm{max}\{0,j-k\}}^{\textrm{min}\{j-1,p\}}
   \binom{p}{t}(-1)^{t}M_{k}(u_j-t)~~.\label{dmnj2}

In DL_POLY_4 the stress tensor due to the reciprocal space, for an
arbitrary :math:`p_\textrm{th}` order multipolar electrostatic
interaction is computed by the formula

.. math::
   :label: virialtensorcomponents_eq

   \begin{aligned}
   V\sigma_{\alpha \beta}^{\textrm{rec}} = \frac{1}{2V_o\epsilon_0\epsilon} \displaystyle \sum_{\mathbf{k}\ne 0}
   \frac{\textrm{exp}(-k^2/4\eta^2)}{k^2} &\left\{ \left
   |S(\mathbf{k})\right|^2 \left[\delta_{\alpha \beta}-2 \left(
   \frac{k^2/4\eta^2 + 1}{k^2}\right) k_{\alpha} k_{\beta}\right] \right. \nonumber \\ 
   & \left. +2S(\mathbf{k}) S_i^{\beta}(-\mathbf{k}) \frac{k_{\alpha}}{k_{\beta}} \right\}~~,
   \end{aligned}

where

.. math:: \mathcal{J}_i^{\textbf{ \ell}}(\mathbf{k}) = \mathcal{M}_i^{\textbf{ \ell}}{\partial}_i^{\textbf{ \ell}}{\textrm{e}}^{\imath \mathbf{k}\cdot \mathbf{r_i}}~~,

.. math:: S_i^{\beta}(-\mathbf{k}) = \sum_{\textbf{ \ell}= \mathbf{0}}^{p} \ell_{\beta} \sum_{i=1}^{N} \mathcal{J}_i^{\textbf{ \ell}}(-\mathbf{k})~~,

and :math:`\textbf{ \ell}= (\ell_1,\ell_2,\ell_3)` .

In DL_POLY_4 the routine ewald_spme_mforces_d computes the reciprocal
space interactions explicitly for simulations with multipoles of maximum
order 2. The routine ewald_spme_mforces handles the general version with
multipoles up to order 4.

The DL_POLY_4 subroutines required to calculate the contributions from
the reciprocal space, in addition to the routines used for the point
charges, are:

#. ``bspgen_mpoles``, in ``spme_container`` evaluates
   equation :eq:`dprodmn_eq` or :eq:`dprodmnsimple_eq`
   to compute the B-splines.

#. ``limit_erfr_deriv`` in ``mpoles_container`` which computes the limit of the
   derivatives of the kernel for the self-interaction term.
   ``limit_erfr_deriv`` is called in ``ewald_spme_mforces``.


.. [1]
   Unlike the other elements of the force field, the electrostatic
   forces are NOT specified in the input FIELD file, but by setting
   appropriate directives in the CONTROL file. See
   Section :ref:`control-file`.

.. [2]
   Strictly speaking, the real space sum ranges over all periodic images
   of the simulation cell, but in the DL_POLY_4 implementation, the
   parameters are chosen to restrict the sum to the simulation cell and
   its nearest neighbours, i.e. the *minimum images* of the cell
   contents.