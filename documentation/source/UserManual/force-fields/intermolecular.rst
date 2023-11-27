.. _intermolecular-potentials:

The Intermolecular Potential Functions
======================================

.. index:: 
   single: potential;metal 
   single: potential;Tersoff 
   single: potential;three-body
   single: potential;four-body 
   single: potential;intramolecular

In this section we outline the two-body, metal, Tersoff, three-body and
four-body potential functions in . An important distinction between
these and intramolecular (bond) forces in DL_POLY_4 is that they are
specified by *atom types* rather than atom indices.

.. _vdw:

Short Ranged (van der Waals) Potentials
---------------------------------------

The short ranged pair forces available in DL_POLY_4 are as follows:

#. 12-6 potential: (\ **12-6**)

   .. math::
      :label: 12-6_vdw_eq

      U(r_{ij}) =
      \left(\frac{A}{r_{ij}^{12}}\right)-\left(\frac{B}{r_{ij}^{6}}\right)

#. Lennard-Jones potential: (\ **lj**)

   .. math::
      :label: lj_vdw_eq

      U(r_{ij}) = 4\epsilon\left[\left
      (\frac{\sigma}{r_{ij}}\right)^{12}-\left(\frac{\sigma}{r_{ij}}\right)^{6}\right]

#. Lennard-Jones cohesive potential :cite:`barrat-99a`: (\
   **ljc**)

   .. math::
      :label: lj_cohesive_vdw_eq

      U(r_{ij}) = 4\epsilon\left[\left
      (\frac{\sigma}{r_{ij}}\right)^{12}-c_{ij}\left(\frac{\sigma}{r_{ij}}\right)^{6}\right]

   This potential has an extra constant to tune the attractive part of
   the potential serving to describe the different cohesiveness between
   different fluids and surfaces in engineering flows models.

#. Lennard Jones, generalised by Frenkel et al.
   :cite:`wang2019`: (\ **ljf**)

   .. math::
      :label: lj_frenkel_eq

      U(r_{ij}) = \left\{
      \begin{array}{lr}
      \varepsilon_{AB}\alpha_{AB}\left[ \left(\frac{\sigma_{AB}}{r_{ij}}\right)^2 -1\right]\left[ \left(\frac{r_{AB}^c}{r_{ij}}\right)^2 -1\right]^2 & r \leq r_{AB}^c \\
      0 & r >  r_{AB}^c \\
      \end{array}
      \right.
        \label{eq:ljf}

   with

   .. math:: \alpha_{AB} = \frac{1}{4}\left(\frac{r_{AB}^c}{\sigma_{AB}}\right)^2\left( \frac{3}{\left(\frac{r_{AB}^c}{\sigma_{AB}}\right)^2-1}\right)^3

#. n-m potential (aka Mie) :cite:`mie-03a,clarke-86a`: (\
   **nm**)

   .. math::
      :label: nm_vdw_eq

      U(r_{ij}) = \frac{E_{o}}{(n-m)}\left[m\left
      (\frac{r_{o}}{r_{ij}}\right)^{n}-n\left(\frac{r_{o}}{r_{ij}}\right)^{m}\right]

#. Buckingham potential: (\ **buck**)

   .. math::
      :label: buck_vdw_eq

      U(r_{ij}) =
      A~\exp\left(-\frac{r_{ij}}{\rho}\right)-\frac{C}{r_{ij}^{6}}

#. Born-Huggins-Meyer potential: (\ **bhm**)

   .. math::
      :label: bhm_vdw_eq

      U(r_{ij}) =
      A~\exp[B(\sigma-r_{ij})]-\frac{C}{r_{ij}^{6}}-\frac{D}{r_{ij}^{8}}

#. Hydrogen-bond (12-10) potential: (\ **hbnd**)

   .. math::
      :label: hbond_vdw_eq

      U(r_{ij}) =
      \left(\frac{A}{r_{ij}^{12}}\right)-\left(\frac{B}{r_{ij}^{10}}\right)

#. Shifted force n-m potential (aka Mie)
   :cite:`mie-03a,clarke-86a`: (\ **snm**)

   .. math::
      :label: shifted_nm_vdw_eq

      \begin{aligned}
      U(r_{ij})=&\frac{\alpha E_{o}}{(n-m)}\left [
      m\beta^{n}\left \{ \left (\frac{r_{o}}{r_{ij}}\right )^{n}-
      \left(\frac{1}{\gamma}\right)^{n}\right \}-
      n\beta^{m}\left \{ \left (\frac{r_{o}}{r_{ij}}\right )^{m}-
      \left(\frac{1}{\gamma}\right)^{m}\right \} \right ]~+~\phantom{xxxx} \nonumber \\
      & \frac{nm\alpha E_{o}}{(n-m)} \left ( \frac{r_{ij}-\gamma r_{o}}{\gamma r_{o}}
      \right )\left\{\left(\frac{\beta}{\gamma}\right
      )^{n}-\left(\frac{\beta}{\gamma}\right )^{m}\right \}~~,\end{aligned}

   with

   .. math::

      \begin{aligned}
      \alpha=&\frac{(n-m)}{[n\beta^{m}(1+(m/\gamma-m-1)/\gamma^{m})-
      m\beta^{n}(1+(n/\gamma-n-1)/\gamma^{n})]} \nonumber \\
      \beta =& \gamma\left( \frac{\gamma^{m+1}-1}{\gamma^{n+1}-1}
      \right)^{\frac{1}{n-m}} \\
      \gamma =& \frac{r_{\rm cut}}{r_{o}}~~. \nonumber\end{aligned}

   This peculiar form has the advantage over the standard shifted n-m
   potential in that both :math:`E_{o}` and :math:`r_{0}` (well depth
   and location of minimum) retain their original values after the
   shifting process.

#. Morse potential: (\ **mors**)

   .. math:: U(r_{ij}) = E_{o}~[\{1-\exp(-k(r_{ij}-r_{o}))\}^{2}-1]
      :label: morse_vdw_eq

#. Shifted Weeks-Chandler-Andersen (WCA) potential
   :cite:`weeks-71a`: (\ **wca**)

   .. math::
      :label: wca_eq

      U(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
      4\epsilon\left[\left(\frac{\sigma}{r_{ij}-\Delta}\right)^{12}-\left(\frac{\sigma}{r_{ij}-\Delta}\right)^{6}\right]
      +\epsilon & r_{ij} < 2^{1 \over 6}~\sigma + \Delta \\
      0 & r_{ij} \ge 2^{1 \over 6}~\sigma + \Delta \end{array} \right. \label{wca}

   The WCA potential is the Lennard-Jones potential truncated at the
   position of the minimum and shifted to eliminate discontinuity
   (includes the effect of excluded volume). It is usually used in
   combination with the FENE, equation :eq:`FENE_bond_eq`, bond
   potential. This implementation allows for a radius shift of up to
   half a :math:`\sigma` (:math:`|\Delta| \le 0.5~\sigma`) with a
   default of zero (:math:`\Delta_{default} = 0`).

#. Standard DPD potential: (\ **dpd**)

   .. math::
      :label: dpd_vdw_eq

      U(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
      \frac{A}{2}~r_{c}~\left(1-\frac{r_{ij}}{r_{c}}\right)^{2} & r_{ij} < r_{c} \\
      0 & r_{ij} \ge r_{c} \end{array} \right.

   It takes the Groot-Warren :cite:`groot-97a` form giving a
   soft and purely repulsive interaction.

#. :math:`n`\ DPD potential: (\ **ndpd**)

   .. math:: 
      
      U(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
      \frac{Ab}{n+1}~r_{c}~\left(1-\frac{r_{ij}}{r_{c}}\right)^{n+1} - \frac{A}{2}~r_{c}~\left(1-\frac{r_{ij}}{r_{c}}\right)^{2} & r_{ij} < r_{c} \\
      0 & r_{ij} \ge r_{c} \end{array} \right.

   It is a modification of the 'standard DPD' Groot-Warren form :cite:`sokhan-23a`, providing 
   additional attraction and more control over repulsion with power index :math:`n` while still 
   remaining soft.

#. 14-7 pair potential :cite:`ponder-10a`: (\ **14-7**)

   .. math:: U(r_{ij}) = \epsilon\left(\frac{1.07}{(r_{ij}/r_{o})+0.07}\right)^{7}\left(\frac{1.12}{(r_{ij}/r_{o})^{7}+0.12}-2\right)
      :label: 14-7_vdw_eq

#. Morse modified :cite:`pedone2006`: (\ **mstw**)

   .. math:: U\left(r_{ij}\right) = E_{0}~\{\left[1-\exp\left(-k\left(r_{ij}-r_{0}\right)\right)\right]^{2}-1\}+\frac{c}{r_{ij}^{12}}
      :label: morse_mod_vdw_eq

#. Rydberg: ( **ryd**)

   .. math:: U\left(r_{ij}\right) = \left(a+br_{ij}\right)\exp\left(-r_{ij}/\rho\right)
      :label: rydberg_vdw_eq

#. Ziegler-Biersack-Littmark (ZBL): :cite:`ziegler1985` (\
   **zbl**)

   .. math::
      :label: zbl_vdw_eq

      U(r_{ij}) = \frac{Z_{1} Z_{2}e^{2}}{4 \pi \varepsilon_{0} \varepsilon_{r}}
                        \sum_{i=1}^{4} b_{i} \exp\left(-c_{i}r/a\right)~~,

   where

   .. math::

      \begin{split}
            a=&~\frac{0.88534 \cdot a_{B}}{Z_{1}^{0.23}+Z_{2}^{0.23}} \\
            b=&~[0.18175,0.50986,0.28022,0.02817] \\
            c=&~[3.1998,0.94229,0.40290,0.20162] \\
        a_{B}=&~ 0.52917721067~\textrm{\AA}~~. \nonumber
          \end{split}

#. ZBL mixed with Morse, :cite:`trachenko2003`: (\ **zbls**)

   .. math:: U\left(r_{ij}\right) = f\left(r_{ij}\right)U_{ZBL}\left(r_{ij}\right)+\left(1-f\left(r_{ij}\right)\right)U_{morse}\left(r_{ij}\right)~~,
      :label: zbl_morse_vdw_eq

   with :math:`f\left(r\right)` defined by

   .. math::

      f\left(r\right) = \left\{
             \begin{array}{ll}
               1-e^{-\left(r_m-r\right)/\xi}/2 & ~~:~~ r<r_m \\
               e^{-\left(r-r_m\right)/\xi}/2 & ~~:~~ r\geq r_m~~. \\
             \end{array}
             \right.

#. ZBL mixed with Buckingham, :cite:`trachenko2003`: (\
   **zblb**)

   .. math:: U\left(r_{ij}\right) = f\left(r_{ij}\right)U_{ZBL}\left(r_{ij}\right)+\left(1-f\left(r_{ij}\right)\right)U_{buckingham}\left(r_{ij}\right)~~,
      :label: zbl_buck_vdw_eq

   with :math:`f\left(r\right)` defined by

   .. math::

      f\left(r\right) = \left\{
             \begin{array}{ll}
               1-e^{-\left(r_m-r\right)/\xi}/2 & ~~:~~ r<r_m \\
               e^{-\left(r-r_m\right)/\xi}/2 & ~~:~~ r\geq r_m~~. \\
             \end{array}
             \right.

#. Lennard-Jones tapered with Mei-Davenport-Fernando taper (MDF),
   :cite:`Raiteri2010`: (\ **mlj**)

   .. math:: U\left(r_{ij}\right) = f\left(r_{ij}\right)U_{LJ}\left(r_{ij}\right)~~,
      :label: lj_mdf_vdw_eq

   where

   .. math::

      f(r)=
           \begin{cases}
             1 & ~~:~~ r<r_i \\
         \frac{\left(r_c-r\right)^3\left(10r_i^2-5r_cr_i-15rr_i+r_c^2+3rr_c+6r^2\right)}{\left(r_c-r_i\right)^5} & ~~:~~ r_i \leq r \leq r_c \\
             0 & ~~:~~ r>r_c~~.\\
           \end{cases}

   :math:`r_c` is set to :math:`r_{vdw}` and controled by **rvdw**.

#. Buckingham tapered with MDF: (\ **mbuc**)

   .. math:: U\left(r_{ij}\right) = f\left(r_{ij}\right)U_{Buckingham}\left(r_{ij}\right)
      :label: buck_mdf_vdw_eq

   See **mlj** for more details.

#. 12-6 Lennard-Jones tapered with Mei-Davenport-Fernando taper (MDF): (\
   **m126**)

   .. math:: U\left(r_{ij}\right) = f\left(r_{ij}\right)U_{12-6}\left(r_{ij}\right)
      :label: 12-6_lj_mdf_vdw_eq

   See **mlj** for more details.

#. Tabulation: (\ **tab**). The potential is defined numerically only.

The parameters defining these potentials are supplied to DL_POLY_4 at
run time (see the description of the FIELD file in
Section :ref:`field-file`). Each atom type in the system
is specified by a unique eight-character label defined by the user. The
pair potential is then defined internally by the combination of two atom
labels.

It is worth noting that some potentials are implemented in an extended
form from their original reference specification. Often this is done by
replacing the :math:`r` argument by :math:`r-r_{o}` to define a surface
softness/hardness width/radius.

As well as the numerical parameters defining the potentials, should also
be provided with a cutoff radius, :math:`r_{\rm vdw}`, which sets a
range limit on the computation of the interactions. It is worth noting
that some interaction come with a hard-wired cutoff in their parameter
sets! Thus any provided cutoff radius, :math:`r_{\rm vdw}`, will be
reset if it is not equal or larger that the largest of these all.
Together with the parameters, the cutoff is used by the subroutine
vdw_generate to construct an interpolation array vvdw for the potential
function over the range 0 to :math:`r_{\rm vdw}`. A second array gvdw is
also calculated, which is related to the potential via the formula:

.. math:: G(r_{ij}) = -r_{ij}\frac{\partial}{\partial r_{ij}}U(r_{ij})~~,

and is used in the calculation of the forces. Both arrays are tabulated
in units of energy. The use of interpolation arrays, rather than the
explicit formulae, makes the routines for calculating the potential
energy and atomic forces very general, and enables the use of user
defined pair potential functions. DL_POLY_4 also allows the user to read
in the interpolation arrays directly from a file (implemented in the
vdw_table_read routine) and the TABLE file
(Section :ref:`table-file`). This is particularly useful if
the pair potential function has no simple analytical description (e.g.
spline potentials).

The force on an atom :math:`j` derived from one of these potentials is
formally calculated with the standard formula:

.. math::
   :label: vdwf_eq

   \underline{f}_{j} = -\frac{1}{r_{ij}}\left[\frac{\partial}{\partial
   r_{ij}}U(r_{ij})\right]\underline{r}_{ij}~~,

where :math:`\underline{r}_{ij} = \underline{r}_{j}-\underline{r}_{i}` . The force on atom
:math:`i` is the negative of this.

The contribution to be added to the atomic virial (for each pair
interaction) is

.. math:: {\cal W} = -\underline{r}_{ij} \cdot \underline{f}_{j}~~.

The contribution to be added to the atomic :index:`stress tensor` is given by

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha}f_{j}^{\beta}~~,
   :label: ast_vdw_eq

where :math:`\alpha` and :math:`\beta` indicate the :math:`x,y,z`
components. The atomic :index:`stress` tensor derived from the pair forces is
symmetric.

Since the calculation of pair potentials assumes a spherical cutoff
(:math:`r_{\rm vdw}`) it is necessary to apply :index:`a<long-ranged corrections;van der Waals>` 
*long-ranged
correction* to the system potential energy and virial. Explicit formulae
are needed for each case and are derived as follows. For two atom types
:math:`a` and :math:`b`, the correction for the potential energy is
calculated via the integral

.. math::

   U_{corr}^{ab} = 2\pi
   \frac{N_{a}N_{b}}{V}\int_{r_{\rm vdw}}^{\infty}g_{ab}(r)U_{ab}(r)r^{2}dr~~,

where :math:`N_{a},N_{b}` are the numbers of atoms of types :math:`a`
and :math:`b` in the system, :math:`V` is the system volume and
:math:`g_{ab}(r)` and :math:`U_{ab}(r)` are the appropriate pair
correlation function and pair potential respectively. It is usual to
assume :math:`g_{ab}(r)=1` for :math:`r>r_{\rm vdw}` . DL_POLY_4
sometimes makes the additional assumption that the repulsive part of the
short ranged potential is negligible beyond :math:`r_{\rm vdw}` .

The correction for the system virial is

.. math::

   {\cal W}_{corr}^{ab} = -2\pi \frac{N_{a}N_{b}}{V}
   \int_{r_{\rm vdw}}^{\infty}g_{ab}(r) \frac{\partial}{\partial
   r}U_{ab}(r)r^{3}dr~~,

where the same approximations are applied.

.. note::
   
   These formulae are based on the assumption that the system
   is reasonably isotropic beyond the cutoff. It is worth noting that the
   14-7 pair potential’s corrections to system energy and virial are solved
   numerically.

In DL_POLY_4 the short ranged forces are calculated by the subroutine
``vdw_forces``. The long-ranged corrections are calculated by routine
``vdw_lrc``. The calculation makes use of the :index:`Verlet<algorithm;Verlet>` 
neighbour list (see above).

Notes on mixing rules for short-ranged interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DL_POLY_4 allows a short cut for mixing some of the explicitly specified
pair interactions for single species of the same type so that
cross-species interactions are generated if unspecified. This is only
possible for the **12-6, lj, dpd, 14-7, wca & ljc** types. The mixing is
derived from the Lennard-Jones style characteristic paramteres for
energy (:math:`\epsilon`) and distance (:math:`\sigma` or :math:`r_{0}`)
terms. The available types of mixing within DL_POLY_4 are borrowed from
:cite:`al-matar-04a`. The rules’ names and formulae are as
follows:

#. Lorentz-Berthelot

   .. math:: \epsilon_{ij} = \sqrt{\epsilon_{i}~\epsilon_{j}}~~;~~\sigma_{ij} = \frac{\sigma_{i}+\sigma_{j}}{2}

#. Fender-Halsey

   .. math:: \epsilon_{ij} = 2 \frac{\epsilon_{i}~\epsilon_{j}}{\epsilon_{i}+\epsilon_{j}}~~;~~\sigma_{ij} = \frac{\sigma_{i}+\sigma_{j}}{2}

#. Hogervorst (good hope)

   .. math:: \epsilon_{ij} = \sqrt{\epsilon_{i}~\epsilon_{j}}~~;~~\sigma_{ij} = \sqrt{\sigma_{i}~\sigma_{j}}

#. Halgren HHG

   .. math:: \epsilon_{ij} = 4 \frac{\epsilon_{i}~\epsilon_{j}}{\left(\epsilon_{i}^{1/2}+\epsilon_{j}^{1/2}\right)^{2}}~~;~~\sigma_{ij} = \frac{\sigma_{i}^{3}+\sigma_{j}^{3}}{\sigma_{i}^{2}+\sigma_{j}^{2}}

#. Waldman-Hagler

   .. math:: \epsilon_{ij} = 2 \sqrt{\epsilon_{i}~\epsilon_{j}} \frac{(\sigma_{i}~\sigma_{j})^{3}}{\sigma_{i}^{6}+\sigma_{j}^{6}}~~;~~\sigma_{ij} = \left(\frac{\sigma_{i}^{6}+\sigma_{j}^{6}}{2}\right)^{\frac{1}{6}}

#. Tang-Toennies

   .. math:: \epsilon_{ij} \sigma_{ij}^{6} = \sqrt{\epsilon_{i} \sigma_{i}^{6}~\epsilon_{j} \sigma_{j}^{6}}~~;~~\epsilon_{ij} \sigma_{ij}^{12}=\left[\frac{\left(\epsilon_{i} \sigma_{i}^{12}\right)^{13}+\left(\epsilon_{j} \sigma_{j}^{12}\right)^{13}}{2}\right]^{13}

#. Functional

   .. math::

      \epsilon_{ij} = \frac{3~\sqrt{\epsilon_{i}~\epsilon_{j}}~(\sigma_{i}~\sigma_{j})^{3}}
      {\sum\limits_{L=0}^{2}{\left[\frac{\left(\sigma_{i}^{3}+\sigma_{j}^{3}\right)^{2}}{4~(\sigma_{i}~\sigma_{i})^{L}}\right]^{\frac{6}{6-2L}}}}~~;~~\sigma_{ij}=\frac{1}{3}\sum\limits_{L=0}^{2}{\left[\frac{\left(\sigma_{i}^{3}+\sigma_{j}^{3}\right)^{2}}{4~(\sigma_{i}~\sigma_{i})^{L}}\right]^{\frac{1}{6-2L}}}

It is woth noting that the :math:`i` and :math:`j` symbols in the
equations for mixing denote atom types (species) and the indices for the
same species interaction parameters are contracted to a single species
index for simplicity.

.. _metal:

Metal Potentials
----------------

.. index:: single: potential;metal

The metal potentials in DL_POLY_4 follow two similar but distinct
formalisms. The first of these is the embedded atom model (EAM)
:cite:`baskes-84a,baskes-86a` and the second is the
Finnis-Sinclair model (FS) :cite:`finnis-84a`. Both are
density dependent potentials derived from density functional theory
(DFT) and describe the bonding of a metal atom ultimately in terms of
the local electronic density. They are suitable for calculating the
properties of metals and metal alloys. The extended EAM (EEAM)
:cite:`hepburn-08a,lau-07a` is a generalisation of the EAM
formalism which can include both EAM and FS type of mixing rules (see
below).

It is worth noting that the same formalism applies to the many-body
perturbation component of the actinide oxide potentials as in
:cite:`cooper-14a`. Thus their many-body component
description is included in this Section.

For single component metals the two main approaches, FS and EAM, are the
same. **However**, they are subtly different in the way they are
extended to handle alloys (see below). It follows that EAM and FS class
potentials cannot be mixed in a single simulation. Furthermore, even for
FS class potentials possessing different analytical forms there is no
agreed procedure for mixing the parameters. Mixing EAM and EEAM
potentials is only possible if the EAM ones are generalised to EEAM form
(see below). The user is, therefore, strongly advised to be consistent
in the choice of potential when modelling alloys.

The general form of the EAM and FS types of potentials is
:cite:`friedel-52a`

.. math::
   :label: um_eq

   U_{metal} = {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} V_{ij}(r_{ij}) +
   \sum_{i=1}^{N} F(\rho_{i})~~, \label{um}

where :math:`F(\rho_{i})` is a functional describing the energy of
embedding an atom in the bulk density, :math:`\rho_{i}`, which is
defined as

.. math:: \rho_{i} = \sum_{j=1, j \ne i}^{N} \rho_{ij}(r_{ij})~~. \label{umd}
   :label: umd_eq

It should be noted that the density is determined by the coordination
number of the atom defined by *pairs* of atoms. This makes the metal
potential dependent on the local density (environmental).
:math:`V_{ij}(r_{ij})` is a pair potential incorporating repulsive
electrostatic and overlap interactions. :math:`N` is the number of
interacting particles in the MD box.

In DL_POLY_4 EAM and thus EEAM can be further generalised to include
two-band (2B) densities :cite:`ackland-03a,ollson-05a`, for
:math:`s`- and :math:`d`-bands,

.. math:: F(\rho_{i})=F^{s}(\rho^{s}_{i})+F^{d}(\rho^{d}_{i})~~, \label{2b}
   :label: 2b_eq

where

.. math:: \rho^{q}_{i} = \sum_{j=1, j \ne i}^{N} \rho^{q}_{ij}(r_{ij})~,~~q=s,d~~, \label{2umd}
   :label: 2umd_eq

instead of just the one, :math:`s`, as in equations :eq:`um_eq`
and :eq:`umd_eq`. These will be referred in the following text as
2BEAM and 2BEEAM. Mixing 2BEAM and EAM and alternatively 2BEEAM and EEAM
potentials is only possible if the single band ones are generalised to
2B forms. The user is, again, reminded to be consistent in the choice of
potential when modelling alloys.

The types of metal potentials available in DL_POLY_4 are as follows:

#. EAM :index:`potential<potential;EAM>`: (\ **eam**) There are no explicit mathematical
   expressions for EAM potentials, so this potential type is read
   exclusively in the form of interpolation arrays from the TABEAM table
   file (as implemented in the metal_table_read routine -
   Section :ref:`tabeam-file`.) The rules for combining
   the potentials from different metals to handle alloys are different
   from the FS class of potentials (see below).

#. EEAM :index:`potential<potential;EEAM>` (\ **eeam**) Similar to EAM above, it is given in the
   form of interpolation arrays from the TABEAM file, but the rules for
   combining the potentials from different metals are different from
   both EAM and FS classes (see below).

#. 2BEAM :index:`potential<potential;2BEAM>` (\ **2beam**) Similar to EEAM for the :math:`s`
   density terms and to EAM for the :math:`d` ones. It is and given in
   the form of interpolation arrays from the TABEAM file, but the rules
   for combining the potentials from different metals are different from
   both EAM, EEAM and FS classes (see below).

#. 2BEEAM :index:`potential<potential;2BEEAM>` (\ **2beeam**) Similar to EEAM for both :math:`s` and
   :math:`d` density terms. It is and given in the form of interpolation
   arrays from the TABEAM file, but the rules for combining the
   potentials from different metals are different from both EAM, EEAM,
   2BEAM and FS classes (see below).

#. Finnis-Sinclair potential :cite:`finnis-84a`: (\ **fnsc**)
   Finnis-Sinclair potential is explicitly analytical. It has the
   following form:

   .. math::

      \begin{aligned}
      V_{ij}(r_{ij})=& \left\{ \begin{array} {l@{\quad:\quad}l}
      (r_{ij}-c)^{2} (c_{0}+c_{1}r_{ij}+c_{2}r_{ij}^{2}) & r_{ij} < c \\
      0 & r_{ij} > c
      \end{array} \right. \nonumber \\
      \rho_{ij}(r_{ij}) =& \left\{ \begin{array} {l@{\quad:\quad}l}
      (r_{ij}-d)^{2} + \beta \displaystyle \frac{(r_{ij}-d)^{3}}{d} & r_{ij} < d \\
      0 & r_{ij} > d
      \end{array} \right. \\
      F(\rho_{i}) =& -A \sqrt{\rho_{i}}~~, \nonumber\end{aligned}

   with parameters: :math:`c_{0}`, :math:`c_{1}`, :math:`c_{2}`,
   :math:`c`, :math:`A`, :math:`d`, :math:`\beta`, both :math:`c` and
   :math:`d` are cutoffs. Since first being proposed a number of
   alternative analytical forms have been proposed, some of which are
   described below. The rules for combining different metal potentials
   to model alloys are different from the EAM potentials (see below).

#. Extended Finnis-Sinclair potential :cite:`dai-06a`: (\
   **exfs**) It has the following form:

   .. math::

      \begin{aligned}
      V_{ij}(r_{ij}) =& \left\{ \begin{array} {l@{\quad:\quad}l}
      (r_{ij}-c)^{2} (c_{0}+c_{1}r_{ij}+c_{2}r_{ij}^{2}+c_{3}r_{ij}^{3}+c_{4}r_{ij}^{4}) & r_{ij} < c \\
      0 & r_{ij} > c
      \end{array} \right. \nonumber \\
      \rho_{ij}(r_{ij}) =& \left\{ \begin{array} {l@{\quad:\quad}l}
      (r_{ij}-d)^{2} + B^{2} (r_{ij}-d)^{4} & r_{ij} < d \\
      0 & r_{ij} > d
      \end{array} \right. \\
      F(\rho_{i}) =& -A \sqrt{\rho_{i}}~~, \nonumber\end{aligned}

   with parameters: :math:`c_{0}`, :math:`c_{1}`, :math:`c_{2}`,
   :math:`c_{3}`, :math:`c_{4}`, :math:`c`, :math:`A`, :math:`d`,
   :math:`B`, both :math:`c` and :math:`d` are cutoffs.

#. Sutton-Chen potential
   :cite:`sutton-90a,sutton-91a,todd-93a`: (\ **stch**) The
   Sutton Chen potential is an analytical potential in the FS class. It
   has the form:

   .. math::

      \begin{aligned}
      V_{ij}(r_{ij}) =& \epsilon \left( \frac{a}{r_{ij}} \right)^{n} \nonumber \\
      \rho_{ij}(r_{ij}) =& \left( \frac{a}{r_{ij}} \right)^{m} \\
      F(\rho_{i}) =& -c \epsilon \sqrt{\rho_{i}}~~, \nonumber\end{aligned}

   with parameters: :math:`\epsilon`, :math:`a`, :math:`n`, :math:`m`,
   :math:`c`. **Note** that the parameter :math:`c` for the mixed
   potential in multi-component allys is irrelevant as outlined in
   :cite:`sutton-91a`!

#. Gupta potential :cite:`cleri-93a`: (\ **gupt**) The Gupta
   potential is another analytical potential in the FS class. It has the
   form:

   .. math::

      \begin{aligned}
      V_{ij}(r_{ij}) =& 2 A \exp \left(-p \frac{r_{ij}-r_{0}}{r_{0}}\right) \nonumber \\
      \rho_{ij}(r_{ij}) =& \exp \left(-2 q_{ij} \frac{r_{ij}-r_{0}}{r_{0}}\right) \\
      F(\rho_{i}) =& -B \sqrt{\rho_{i}}~~, \nonumber\end{aligned}

   with parameters: :math:`A`, :math:`r_{0}`, :math:`p`, :math:`B`,
   :math:`q_{ij}`.

#. Many body perturbation component potential
   :cite:`cooper-14a`: (\ **mbpc**) This component is another
   analytical potential in the FS class which two body part may be
   defined by a matching van der Waals potential in the vdw section of
   the FIELD file. It has the form:

   .. math::

      \begin{aligned}
      V_{ij}(r_{ij}) =& 0 \nonumber \\
      \rho_{ij}(r_{ij})=& \left(\frac{a}{r_{ij}^{m}}\right) \frac{1}{2}\left[1+{\rm erf}\left(\alpha(r_{ij}-r_{\rm o})\right)\right] \\
      F(\rho_{i}) =& -\epsilon \sqrt{\rho_{i}}~~, \nonumber\end{aligned}

   with parameters: :math:`\epsilon`, :math:`a`, :math:`m`,
   :math:`\alpha` and :math:`r_{\rm o}`.

   .. note::
      
      The parameters :math:`\alpha` and :math:`r_{\rm o}`
      must be the same for all defined potentials of this type. DL_POLY_4
      will set :math:`\alpha={\rm Max}(0,\alpha_{pq})` and
      :math:`r_{\rm o}={\rm Max}(0,r_{{\rm o}\_pq})` for all defined
      interactions of this type between species :math:`p` and :math:`q`. If
      after this any is left undefined, i.e. zero, the undefined entities
      will be set to their defaults: :math:`\alpha=20` and
      :math:`r_{\rm o}={\rm Min}(1.5,0.2~r_{\rm cut})`.

All of these metal potentials can be decomposed into pair contributions
and thus fit within the general tabulation scheme of , where they are
treated as pair interactions (though note that the metal cutoff,
:math:`r_{\rm met}` has nothing to do with short ranged cutoff,
:math:`r_{\rm vdw}`). DL_POLY_4 calculates this potential in two stages:
the first calculates the local density, :math:`\rho_{i}`, for each atom;
and the second calculates the potential energy and forces. Interpolation
arrays, vmet, gmet and fmet (metal_generate, metal_table_read) are used
in both these stages in the same spirit as in the van der Waals
interaction calculations.

The total force :math:`\underline{f}_{k}^{tot}` on an atom :math:`k` derived
from this potential is calculated in the standard way:

.. math:: \underline{f}_{k}^{tot} = -\underline{\nabla}_{k} U_{metal}~~.

We rewrite the EAM/FS potential, equation :eq:`um_eq`, as

.. math::

   \begin{aligned}
   U_{metal} =& U_{1} + U_{2} \nonumber \\
   U_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} V_{ij}(r_{ij}) \\
   U_{2} =& \sum_{i=1}^{N} F(\rho_{i})~~, \nonumber
   \end{aligned}

where :math:`\underline{r}_{ij} = \underline{r}_{j}-\underline{r}_{i}` . The force on atom
:math:`k` is the sum of the derivatives of :math:`U_{1}` and
:math:`U_{2}` with respect to :math:`\underline{r_{k}}`, which is recognisable
as a sum of pair forces:

.. math::

   \begin{aligned}
   -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& -{1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N}
   \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \underline{r_{k}}} =
   \sum_{j=1,j \ne k}^{N} \frac{\partial V_{kj}(r_{kj})}{\partial r_{kj}} \frac{\underline{r_{kj}}}{r_{kj}} \nonumber \\
   -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& -\sum_{i=1}^{N} \frac{\partial F}{\partial \rho_{i}}
   \sum_{j \ne i}^{N} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \underline{r_{k}}} \\
   =& -\sum_{i=1,i \ne k}^{N} \frac{\partial F}{\partial \rho_{i}} \frac{\partial \rho_{ik}(r_{ik})}{\partial r_{ik}}
   \frac{\partial r_{ik}}{\partial \underline{r_{k}}} - \sum_{j=1,j \ne k}^{N} \frac{\partial F}{\partial \rho_{k}}
   \frac{\partial \rho_{kj}(r_{kj})}{\partial r_{kj}} \frac{\partial r_{kj}}{\partial \underline{r_{k}}} \nonumber \\
   =& \sum_{j=1,j \ne k}^{N} \left( \frac{\partial F}{\partial \rho_{k}} + \frac{\partial F}{\partial \rho_{j}} \right)
   \frac{\partial \rho_{kj}(r_{kj})}{\partial r_{kj}} \frac{\underline{r_{kj}}}{r_{kj}}~~. \nonumber
   \end{aligned}

#. | EAM force
   | The same as shown above. However, it is worth noting that the
     generation of the force arrays from tabulated data (implemented in
     the metal_table_derivatives routine) is done using a five point
     interpolation procedure.

#. | EEAM force
   | Information the same as that for EAM.

#. | 2BEAM force
   | Information the same as that for EAM. However, as there is a second
     embedding contribution from the extra band complexity:
     :math:`U_{2}=U^{s}_{2}+U^{d}_{2}` !

#. | 2BEEAM force
   | Information the same as that for EAM. However, as there is a second
     embedding contribution from the extra band complexity:
     :math:`U_{2}=U^{s}_{2}+U^{d}_{2}` !

#. Finnis-Sinclair force

   .. math::

      \begin{aligned}
      -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& \sum_{j=1,j \ne k}^{N} \left\{
      2 (r_{kj}-c) (c_{0}+c_{1}r_{kj}+c_{2}r_{kj}^{2}) +
      (r_{kj}-c)^{2} (c_{1}+2c_{2}r_{kj}) \right\} \frac{\underline{r_{kj}}}{r_{kj}} \nonumber \\
      -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& -\sum_{j=1,j \ne k}^{N}
      {A \over 2} \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \left\{ 2(r_{kj}-d) + 3 \beta \frac{(r_{kj}-d)^{2}}{d} \right\} \frac{\underline{r_{kj}}}{r_{kj}}~~.\end{aligned}

#. Extended Finnis-Sinclair force

   .. math::

      \begin{aligned}
      -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& \sum_{j=1,j \ne k}^{N} \left\{
      2 (r_{kj}-c) (c_{0}+c_{1}r_{kj}+c_{2}r_{kj}^{2}+c_{3}r_{kj}^{3}+c_{4}r_{kj}^{4}) + \right. \nonumber \\
       & \phantom{xxxxxxx} \left. (r_{kj}-c)^{2} (c_{1}+2c_{2}r_{kj}+3c_{3}r_{kj}^{2}+4c_{4}r_{kj}^{3}) \right\}
      \frac{\underline{r_{kj}}}{r_{kj}} \\
      -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& -\sum_{j=1,j \ne k}^{N}
      {A \over 2} \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \left\{ 2(r_{kj}-d) + 4 B^{2} (r_{kj}-d)^{3}\right\} \frac{\underline{r_{kj}}}{r_{kj}}~~. \nonumber\end{aligned}

#. Sutton-Chen force

   .. math::

      \begin{aligned}
      -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& -\sum_{j=1,j \ne k}^{N} n \epsilon
      \left( \frac{a}{r_{kj}} \right)^{n} \frac{\underline{r_{kj}}}{r_{kj}} \nonumber \\
      -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& \sum_{j=1,j \ne k}^{N} \frac{m c \epsilon}{2}
      \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \left( \frac{a}{r_{kj}} \right)^{m} \frac{\underline{r_{kj}}}{r_{kj}}~~.\end{aligned}

#. Gupta force

   .. math::

      \begin{aligned}
      -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& -\sum_{j=1,j \ne k}^{N} \frac{2 A p}{r_{0}}
      \exp \left( -p \frac{r_{kj}-r_{0}}{r_{0}} \right) \frac{\underline{r_{kj}}}{r_{kj}} \nonumber \\
      -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& \sum_{j=1,j \ne k}^{N} \frac{B q_{kj}}{r_{0}}
      \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \exp \left( -2 q_{kj} \frac{r_{kj}-r_{0}}{r_{0}} \right) \frac{\underline{r_{kj}}}{r_{kj}}~~.\end{aligned}

#. Many body perturbation component potential force

   .. math::

      \begin{aligned}
      -\frac{\partial U_{1}}{\partial \underline{r_{k}}} =& 0 \nonumber \\
      -\frac{\partial U_{2}}{\partial \underline{r_{k}}} =& \sum_{j=1,j \ne k}^{N} \frac{m \epsilon}{2}
      \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \frac{a}{r_{kj}^{m}} \frac{\underline{r_{kj}}}{r_{kj}}~~.\end{aligned}

With the metal forces thus defined the contribution to be added to the
atomic virial *from each atom pair* is then

.. math:: {\cal W} = -\underline{r}_{ij} \cdot \underline{f}_{j}~~,

which equates to:

.. math::

   \begin{aligned}
   \Psi =& 3 V \frac{\partial U}{\partial V} \nonumber \\
   \Psi =& {3 \over 2} V \sum_{i=1}^{N} \sum_{j \ne i}^{N}
   \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial V} +
   3 V \sum_{i=1}^{N} \frac{\partial F(\rho_{i})}{\partial \rho_{i}} \frac{\partial \rho_{i}}{\partial V}
   = \Psi_{1} + \Psi_{2} \nonumber \\
   & \frac{\partial r_{ij}}{\partial V} = \frac{\partial V^{1/3}s_{ij}}{\partial V} =
   {1 \over 3} V^{-2/3}s_{ij} = \frac{r_{ij}}{3 V} \nonumber \\
   \Psi_{1} &=& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} \\
   & \frac{\partial \rho_{i}}{\partial V} = \frac{\partial }{\partial V} \sum_{j=1, j \ne i}^{N} \rho_{ij}(r_{ij}) =
   \sum_{j=1, j \ne i}^{N} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial V} =
   \frac{1}{3 V} \sum_{j=1, j \ne i}^{N} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} \nonumber \\
   \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} \left( \frac{\partial F(\rho_{i})}{\partial \rho_{i}} +
   \frac{\partial F(\rho_{j})}{\partial \rho_{j}} \right) \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} r_{ij}~~. \nonumber\end{aligned}

#. | EAM virial
   | The same as above.

#. | EEAM virial
   | The same as above.

#. | 2BEAM virial
   | The same as above but with a second embedding contribution from the
     extra band complexity: :math:`\Psi_{2}=\Psi^{s}_{2}+\Psi^{d}_{2}` !

#. | 2BEEAM virial
   | The same as above but with a second embedding contribution from the
     extra band complexity: :math:`\Psi_{2}=\Psi^{s}_{2}+\Psi^{d}_{2}` !

#. Finnis-Sinclair virial

   .. math::

      \begin{aligned}
      \Psi_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N}
      \left\{ 2 (r_{ij}-c) (c_{0}+c_{1}r_{ij}+c_{2}r_{ij}^{2}) +
      (r_{ij}-c)^{2} (c_{1}+2c_{2}r_{ij}) \right\} r_{ij} \nonumber \\
      \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N}
      {A \over 2} \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \left\{ 2(r_{ij}-d) + 3 \beta \frac{(r_{ij}-d)^{2}}{d} \right\} r_{ij}a~~.
      \end{aligned}

#. Extended Finnis-Sinclair virial

   .. math::

      \begin{aligned}
      \Psi_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N}
      \left\{ 2 (r_{ij}-c) (c_{0}+c_{1}r_{ij}+c_{2}r_{ij}^{2}+c_{3}r_{ij}^{3}+c_{4}r_{ij}^{4}) + \right. \nonumber \\
      & \phantom{xxxxxxxx} \left. (r_{ij}-c)^{2} (c_{1}+2c_{2}r_{ij}+3c_{3}r_{ij}^{2}+4c_{4}r_{ij}i^{3}) \right\} r_{ij} \\
      \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N}
      {A \over 2} \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \left\{ 2(r_{ij}-d) + 4 B^{2} (r_{ij}-d)^{3} \right\} r_{ij}a~~. \nonumber
      \end{aligned}

#. Sutton-Chen virial

   .. math::

      \begin{aligned}
      \Psi_{1} =& -{1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} n \epsilon \left( \frac{a}{r_{ij}} \right)^{n} \nonumber \\
      \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} \frac{m c \epsilon}{2} \left( \frac{\partial F(\rho_{i})}{\partial \rho_{i}} +
      \frac{\partial F(\rho_{j})}{\partial \rho_{j}} \right) \left( \frac{a}{r_{ij}} \right)^{m}~~.
      \end{aligned}

#. Gupta virial

   .. math::

      \begin{aligned}
      \Psi_{1} =& -\sum_{i=1}^{N} \sum_{j \ne i}^{N}
      \frac{A p}{r_{0}} \exp \left( -p \frac{r_{ij}-r_{0}}{r_{0}} \right) r_{ij} \nonumber \\
      \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} \frac{B q_{ij}}{r_{0}}
      \left( {1 \over \sqrt{\rho_{k}}} + {1 \over \sqrt{\rho_{j}}} \right)
      \exp \left( -2 q_{ij} \frac{r_{ij}-r_{0}}{r_{0}} \right) r_{ij}~~.
      \end{aligned}

#. Many body perturbation component virial

   .. math::

      \begin{aligned}
      \Psi_{1} =& 0 \nonumber \\
      \Psi_{2} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{N} \frac{m \epsilon}{2} \left( \frac{\partial F(\rho_{i})}{\partial \rho_{i}} +
      \frac{\partial F(\rho_{j})}{\partial \rho_{j}} \right) \frac{a}{r_{ij}^{m}}~~.
      \end{aligned}

The contribution to be added to the atomic :index:`stress tensor` is given by

.. math:: \sigma^{\alpha \beta} = r_{ij}^{\alpha} f_{j}^{\beta}~~,

where :math:`\alpha` and :math:`\beta` indicate the :math:`x,y,z`
components. The atomic stress tensor is symmetric.

.. index:: single: long-ranged corrections;metal
   
The long-ranged correction for the DL_POLY_4 metal potential is in two
parts. Firstly, by analogy with the short ranged potentials, the
correction to the local density is

.. math::

   \begin{aligned}
   \rho_{i} =& \sum_{j=1, j \ne i}^{\infty} \rho_{ij}(r_{ij}) \nonumber \\
   \rho_{i} =& \sum_{j=1, j \ne i}^{r_{ij}<r_{\rm met}} \rho_{ij}(r_{ij}) +
   \sum_{j=1, j \ne i}^{r_{ij} \ge r_{\rm met}} \rho_{ij}(r_{ij}) =
   \rho_{i}^{o} + \delta \rho_{i} \\
   \delta \rho_{i} =& 4 \pi \bar{\rho} \int_{r_{\rm met}}^{\infty} \rho_{ij}(r) dr~~, \nonumber
   \end{aligned}

where :math:`\rho_{i}^{o}` is the uncorrected local density and
:math:`\bar{\rho}` is the *mean particle density*. Evaluating the
integral part of the above equation yields:

#. | EAM density correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | EEAM density correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEAM density correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEAM density correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | Finnis-Sinclair density correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. | Extended Finnis-Sinclair density correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. Sutton-Chen density correction

   .. math::

      \delta \rho_{i} = \frac{4 \pi \bar{\rho} a^{3}}{(m-3)}
      \left( \frac{a}{r_{\rm met}} \right)^{m-3}~~.

#. Gupta density correction

   .. math::

      \delta \rho_{i} = \frac{2 \pi \bar{\rho} r_{0}}{q_{ij}}
      \left[ r_{\rm met}^{2} + 2 r_{\rm met} \left(\frac{r_{0}}{q_{ij}}\right) +
      2 \left(\frac{r_{0}}{q_{ij}}\right)^{2} \right]
      \exp \left( -2 q_{ij} \frac{r_{\rm met}-r_{0}}{r_{0}}\right)~~.

#. Many body perturbation component density correction

   .. math::

      \delta \rho_{i} = \frac{4 \pi \bar{\rho}}{(m-3)}
      \frac{a}{r_{\rm met}^{m-3}}~~.

The density correction is applied immediately after the local density is
calculated. The pair term correction is obtained by analogy with the
short ranged potentials and is

.. math::

   \begin{aligned}
   U_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{\infty} V_{ij}(r_{ij}) \nonumber \\
   U_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{r_{ij}<r_{\rm met}} V_{ij}(r_{ij}) +
   {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{r_{ij} \ge r_{\rm met}} V_{ij}(r_{ij}) =
   U_{1}^{o} + \delta U_{1} \nonumber \\
   \delta U_{1} =& 2 \pi N \bar{\rho} \int_{r_{\rm met}}^{\infty} V_{ij}(r) r^{2} dr \nonumber \\
   U_{2} =& \sum_{i=1}^{N} F(\rho_{i}^{0} + \delta \rho_{i}) \\
   U_{2} =& \sum_{i=1}^{N} F(\rho_{i}^{0}) +
   \sum_{i=1}^{N} \frac{\partial F(\rho_{i})_{0}}{\partial \rho_{i}} \delta \rho_{i} =
   U_{2}^{0} + \delta U_{2} \nonumber \\
   \delta U_{2} =& 4 \pi \bar{\rho} \sum_{i=1}^{N} \frac{\partial F(\rho_{i})_{0}}{\partial \rho_{i}}
   \int_{r_{\rm met}}^{\infty} \rho_{ij}(r) r^{2} dr~~. \nonumber
   \end{aligned}

.. note:: 

   :math:`\delta U{2}` is not required if :math:`\rho_{i}`
   has already been corrected.

Evaluating the integral part of the above
equations yields:

#. | EAM energy correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | EEAM energy correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEAM energy correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEEAM energy correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | Finnis-Sinclair energy correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. | Extended Finnis-Sinclair energy correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. Sutton-Chen energy correction

   .. math::

      \begin{aligned}
      \delta U_{1} =& \frac{2 \pi N \bar{\rho} \epsilon a^{3}}{(n-3)}
      \left( \frac{a}{r_{\rm met}} \right)^{n-3} \nonumber \\
      \delta U_{2} =& -\frac{4 \pi \bar{\rho} a^{3}}{(m-3)} \left( \frac{a}{r_{\rm met}} \right)^{m-3}
      \left< \frac{N c \epsilon}{2\sqrt{\rho_{i}^{0}}} \right>~~.
      \end{aligned}

#. Gupta energy correction

   .. math::

      \begin{aligned}
      \delta U_{1} =& \frac{4 \pi N \bar{\rho} A r_{0}}{p}
      \left[ r_{\rm met}^{2} + 2 r_{\rm met} \left(\frac{r_{0}}{p}\right) +
      2 \left(\frac{r_{0}}{p}\right)^{2} \right] \times \nonumber \\
      & \exp \left( -p \frac{r_{\rm met}-r_{0}}{r_{0}}\right) \nonumber \\
      \delta U_{2} =& -\frac{2 \pi \bar{\rho} r_{0}}{q_{ij}}
      \left[ r_{\rm met}^{2} + 2 r_{\rm met} \left(\frac{r_{0}}{q_{ij}}\right) +
      2 \left(\frac{r_{0}}{q_{ij}}\right)^{2} \right] \times \\
      & \exp \left( -2 q_{ij} \frac{r_{\rm met}-r_{0}}{r_{0}}\right)
      \left< \frac{N B}{2\sqrt{\rho_{i}^{0}}} \right>~~. \nonumber
      \end{aligned}

#. Many body perturbation component energy correction

   .. math::

      \begin{aligned}
      \delta U_{1} =& 0 \nonumber \\
      \delta U_{2} =& -\frac{4 \pi \bar{\rho}}{(m-3)} \frac{a}{r_{\rm met}^{m-3}}
      \left< \frac{N \epsilon}{2\sqrt{\rho_{i}^{0}}} \right>~~.
      \end{aligned}

To estimate the virial correction we assume the corrected local
densities are constants (i.e. independent of distance - at least beyond
the range :math:`r_{\rm met}`). This allows the virial correction to be
computed by the methods used in the short ranged potentials:

.. math::

   \begin{aligned}
   \Psi_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{\infty}
   \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} \nonumber \\
   \Psi_{1} =& {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{r_{ij}<r_{\rm met}}
   \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} +
   {1 \over 2} \sum_{i=1}^{N} \sum_{j \ne i}^{r_{ij} \ge r_{\rm met}}
   \frac{\partial V_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} = \Psi_{1}^{0} + \delta \Psi_{1} \nonumber \\
   \delta \Psi_{1} =& 2 \pi N \bar{\rho} \int_{r_{\rm met}}^{\infty}
   \frac{\partial V_{ij}(r)}{\partial r_{ij}} r^{3} dr \nonumber \\
   \Psi_{2} =& \sum_{i=1}^{N} \frac{\partial F(\rho_{i})}{\partial \rho_{i}}
   \sum_{j \ne i}^{\infty} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} \\
   \Psi_{2} =& \sum_{i=1}^{N} \frac{\partial F(\rho_{i})}{\partial \rho_{i}}
   \sum_{j \ne i}^{r_{ij}<r_{\rm met}} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} +
   \sum_{i=1}^{N} \frac{\partial F(\rho_{i})}{\partial \rho_{i}}
   \sum_{j \ne i}^{r_{ij} \ge r_{\rm met}} \frac{\partial \rho_{ij}(r_{ij})}{\partial r_{ij}} r_{ij} = \Psi_{2}^{0} + \delta \Psi_{2} \nonumber \\
   \delta \Psi_{2} =& 4 \pi \bar{\rho} \sum_{i=1}^{N} \frac{\partial F(\rho_{i})}{\partial \rho_{i}}
   \int_{r_{\rm met}}^{\infty} \frac{\partial \rho_{ij}(r)}{\partial r} r^{3} dr~~. \nonumber
   \end{aligned}

Evaluating the integral part of the above equations yields:

#. | EAM virial correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | EEAM virial correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEAM virial correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | 2BEEAM virial correction
   | No long-ranged corrections apply beyond :math:`r_{\rm met}`.

#. | Finnis-Sinclair virial correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. | Extended Finnis-Sinclair virial correction
   | No long-ranged corrections apply beyond cutoffs :math:`c` and
     :math:`d`.

#. Sutton-Chen virial correction

   .. math::

      \begin{aligned}
      \delta \Psi_{1} =& -n \frac{2 \pi N \bar{\rho} \epsilon a^{3}}{(n-3)}
      \left( \frac{a}{r_{\rm met}} \right)^{n-3} \nonumber \\
      \delta \Psi_{2} =& m \frac{4 \pi \bar{\rho} a^{3}}{(m-3)} \left( \frac{a}{r_{\rm met}} \right)^{m-3}
      \left< \frac{N c \epsilon}{2\sqrt{\rho_{i}^{0}}} \right>~~.\end{aligned}

#. Gupta virial correction

   .. math::

      \begin{aligned}
      \delta \Psi_{1} =& -\frac{p}{r_{0}} \frac{4 \pi N \bar{\rho} A r_{0}}{p}
      \left[ r_{\rm met}^{3} + 3 r_{\rm met}^{2} \left(\frac{r_{0}}{p}\right) +
      6 r_{\rm met} \left(\frac{r_{0}}{p}\right)^{2} + 6 \left(\frac{r_{0}}{p}\right)^{3} \right] \times \nonumber \\
      & \exp \left( -p \frac{r_{\rm met}-r_{0}}{r_{0}}\right) \nonumber \\
      \delta \Psi_{2} =& \frac{q_{ij}}{r_{0}} \frac{2 \pi \bar{\rho} r_{0}}{q_{ij}}
      \left[ r_{\rm met}^{3} + 3 r_{\rm met}^{2} \left(\frac{r_{0}}{q_{ij}}\right) +
      6 r_{\rm met} \left(\frac{r_{0}}{q_{ij}}\right)^{2} + 6 \left(\frac{r_{0}}{q_{ij}}\right)^{3} \right] \times \\
      & \exp \left( -2 q_{ij} \frac{r_{\rm met}-r_{0}}{r_{0}}\right)
      \left< \frac{N B}{2\sqrt{\rho_{i}^{0}}} \right>~~. \nonumber
      \end{aligned}

#. Many body perturbation component virial correction

   .. math::

      \begin{aligned}
      \delta \Psi_{1} =& 0 \nonumber \\
      \delta \Psi_{2} =& m \frac{4 \pi \bar{\rho}}{(m-3)} \frac{a}{r_{\rm met}^{m-3}}
      \left< \frac{N \epsilon}{2\sqrt{\rho_{i}^{0}}} \right>~~.
      \end{aligned}

In the energy and virial corrections we have used the approximation:

.. math:: \sum_{i}^{N}\rho_{i}^{-1/2} = \frac{N}{<\rho_{i}^{1/2}>}~~,

where :math:`<\rho_{i}^{1/2}>` is regarded as a constant of the system.

In DL_POLY_4 the metal forces are handled by the routine ``metal_forces``.
The local density is calculated by the routines ``metal_ld_collect_eam``,
``metal_ld_collect_fst``, ``metal_ld_compute``, ``metal_ld_set_halo`` and
``metal_ld_export``. The long-ranged corrections are calculated by
``metal_lrc``. Reading and generation of EAM table data from TABEAM is
handled by ``metal_table_read`` and ``metal_table_derivatives``.

Notes on the Treatment of Alloys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The distinction to be made between EAM and FS potentials with regard to
alloys concerns the mixing rules for unlike interactions. Starting with
equations :eq:`um_eq` and :eq:`umd_eq`, it is clear that we
require mixing rules for terms :math:`V_{ij}(r_{ij})` and
:math:`\rho_{ij}(r_{ij})` when atoms :math:`i` and :math:`j` are of
different kinds. Thus two different metals :math:`A` and :math:`B` we
can distinguish 4 possible variants of each:

.. math::

   V^{AA}_{ij}(r_{ij}),~V^{BB}_{ij}(r_{ij}),~V^{AB}_{ij}(r_{ij}),
   ~V^{BA}_{ij}(r_{ij})

and

.. math::

   \rho^{AA}_{ij}(r_{ij}),~\rho^{BB}_{ij}(r_{ij}),~\rho^{AB}_{ij}(r_{ij}),
   ~\rho^{BA}_{ij}(r_{ij})~~.

These forms recognise that the contribution of a type :math:`A` atom to
the potential of a type :math:`B` atom may be different from the
contribution of a type :math:`B` atom to the potential of a type
:math:`A` atom. In both EAM :cite:`johnson-89a` and FS
:cite:`sutton-91a` cases it turns out that

.. math:: V^{BA}_{ij}(r_{ij})=V^{BA}_{ij}(r_{ij})~~,

though the mixing rules are different in each case (\ **beware!**). This
has the following implications to densities of mixtures for different
potential frameworks:

-  EAM case - it is required that :cite:`johnson-89a`:

   .. math::

      \begin{aligned}
      \rho^{AB}_{ij}(r_{ij})=&\rho^{BB}_{ij}(r_{ij}) \nonumber \\
      \rho^{BA}_{ij}(r_{ij})=&\rho^{AA}_{ij}(r_{ij})~~,
      \end{aligned}

   which means that an atom of type :math:`A` contributes the same
   density to the environment of an atom of type :math:`B` as it does to
   an atom of type :math:`A`, and *vice versa*.

-  EEAM case - all densities can be different
   :cite:`hepburn-08a,lau-07a`:

   .. math::

      \rho^{AA}_{ij}(r_{ij}) \neq \rho^{BB}_{ij}(r_{ij}) \neq \rho^{AB}_{ij}(r_{ij}) \neq
      \rho^{BA}_{ij}(r_{ij})~~!

-  2BEAM case - similarly to the EAM case it is required that:

   .. math::

      \begin{aligned}
      {\rho^{d}_{ij}}^{AB}(r_{ij})=&{\rho^{d}_{ij}}^{BB}(r_{ij}) \nonumber \\
      {\rho^{d}_{ij}}^{BA}(r_{ij})=&{\rho^{d}_{ij}}^{AA}(r_{ij})~~,\end{aligned}

   for the :math:`d`-band densities, whereas for the :math:`s`-band
   ones:

   .. math:: {\rho^{s}_{ij}}^{BA}(r_{ij}) = {\rho^{s}_{ij}}^{AB}(r_{ij})~~,

   which means that an atom of type :math:`A` contributes the same
   :math:`s` density to the environment of an atom of type :math:`B` as
   an atom of type :math:`B` to an environment of an atom of type
   :math:`A`. However, in general:

   .. math:: {\rho^{s}_{ij}}^{AA}(r_{ij}) \neq {\rho^{s}_{ij}}^{BB}(r_{ij}) \neq {\rho^{s}_{ij}}^{AB}(r_{ij})~~. \\

-  2BEEAM case - similarly to the EEAM case all :math:`s` and :math:`d`
   densities can be different:

   .. math::

      \begin{aligned}
      {\rho^{s}_{ij}}^{AA}(r_{ij}) \neq {\rho^{s}_{ij}}^{BB}(r_{ij}) \neq
      {\rho^{s}_{ij}}^{AB}(r_{ij}) \neq {\rho^{s}_{ij}}^{BA}(r_{ij}) \nonumber \\
      {\rho^{d}_{ij}}^{AA}(r_{ij}) \neq {\rho^{d}_{ij}}^{BB}(r_{ij}) \neq
      {\rho^{d}_{ij}}^{AB}(r_{ij}) \neq {\rho^{d}_{ij}}^{BA}(r_{ij}) ~~.\end{aligned}

-  FS case - here a different rule applies
   :cite:`sutton-91a`:

   .. math:: \rho^{AB}_{ij}(r_{ij})=(\rho^{AA}_{ij}(r_{ij})~\rho^{BB}_{ij}(r_{ij}))^{1/2}

   so that atoms of type :math:`A` and :math:`B` contribute the same
   densities to each other, but not to atoms of the same type.

The above rules have the following consequences to the specifications of
these potentials in the DL_POLY_4 FIELD file for an alloy composed of
:math:`n` different metal atom types both the EAM types and FS types of
potentials require the specification of :math:`n(n+1)/2` pair functions
:math:`V^{AB}_{ij}(r_{ij})`. However, the its only the simple EAM type
together with all the FS types that require only :math:`n` density
functions :math:`\rho^{AA}_{ij}(r_{ij})`, whereas the EEAM class
requires all the cross functions :math:`\rho^{AB}_{ij}(r_{ij})` possible
or :math:`n^{2}` in total! In addition to the :math:`n(n+1)/2` pair
functions and :math:`n` or :math:`n^{2}` density functions both the EAM
and EEAM potentials require further specification of :math:`n`
functional forms of the density dependence (i.e. the embedding function
:math:`F(\rho_{i})` in equation :eq:`um_eq`. The matter is
further complicated when the 2BEAM type of potential is used with the
extra specification of :math:`n` embedding functions and
:math:`n(n+1)/2` density functions for the :math:`s`-band. Similarly, in
the 2BEEAM an extra :math:`n` embedding functions and :math:`n^{2}`
density functions for the :math:`s`-band are required.

It is worth noting that in the 2BEAM and 2BEEAM the :math:`s`-band
contribution is usually only for the alloy component, so that local
concentrations of a single element revert to the standard EAM or EEAM!
In such case, the densities functions must be zeroed in the DL_POLY_4
TABEAM file.

For EAM, EEAM, 2BEAM and 2BEEAM potentials all the functions are
supplied in tabular form via the table file TABEAM (see
section :ref:`tabeam-file`) to which DL_POLY_4 is
redirected by the FIELD file data. The FS potentials are defined via the
necessary parameters in the FIELD file.

.. _tersoff:

Tersoff Potentials
------------------

.. index:: single: potential;Tersoff

The Tersoff :cite:`tersoff-89a` potential is is a bond-order
potential, developed to be used in multi-component covalent systems by
an effective coupling of two-body and higher many-body correlations into
one model. The central idea is that in real systems, the strength of
each bond depends on the local environment, i.e. an atom with many
neighbors forms weaker bonds than an atom with few neighbors.
Effectively, it is a pair potential the strength of which depends on the
environment. At the present there are two versions of this potential
available in : **ters** and **kihs**. In these particular
implementations **ters** has 11 atomic and 2 bi-atomic parameters
whereas **kihs** :cite:`kumagai-07a` has 16 atomic
parameters. The energy is modelled as a sum of pair-like interactions,
where the coefficient of the attractive term in the pair-like potential
(which plays the role of a bond order) depends on the local environment
giving a many-body potential.

The form of the Tersoff potential is: (\ **ters**)

.. math:: U_{ij} = f_{C}(r_{ij})~[f_{R}(r_{ij}) - \gamma_{ij}~f_{A}(r_{ij})]~~,

where :math:`f_{R}` and :math:`f_{A}` are the repulsive and attractive
pair potential respectively:

.. math::

   f_{R}(r_{ij}) = A_{ij}~\exp(- a_{ij}~r_{ij})~~,~~
   f_{A}(r_{ij}) = B_{ij}~\exp(- b_{ij}~r_{ij})

and :math:`f_{C}` is a smooth cutoff function with parameters :math:`R`
and :math:`S` so chosen that to include the first-neighbor shell:

-  **ters**:

   .. math::

      f_{C}(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
      1 & r_{ij} < R_{ij} \\
      \frac{1}{2} + \frac{1}{2} \cos \left[\pi~\frac{r_{ij}-R_{ij}}{S_{ij}-R_{ij}}\right]
      & R_{ij} < r_{ij} < S_{ij} \\
      0 & r_{ij} > S_{ij}
      \end{array} \right.

-  **kihs** - here :math:`f_{C}` is modified to a have continuous
   second-order differential:

   .. math::

      f_{C}(r_{ij}) = \left\{ \begin{array} {l@{\quad:\quad}l}
      1 & r_{ij} < R_{ij} \\
      \frac{1}{2} + \frac{9}{16} \cos \left[\pi~\frac{r_{ij}-R_{ij}}{S_{ij}-R_{ij}}\right]
      - \frac{1}{16} \cos \left[3 \pi~\frac{r_{ij}-R_{ij}}{S_{ij}-R_{ij}}\right]
      & R_{ij} < r_{ij} < S_{ij} \\
      0 & r_{ij} > S_{ij}~~.
      \end{array} \right.

:math:`\gamma_{ij}` expresses a dependence that can accentuate or
diminish the attractive force relative to the repulsive force, according
to the local environment, such that:

-  **ters**:

   .. math::

      \begin{aligned}
      \gamma_{ij} &=& \chi_{ij}~(1 + {\beta_{i}}^{\eta_{i}}~{\cal{L}}_{ij}^{\eta_{i}})^{-\frac{1}{2\eta_{i}}} \nonumber \\
      {\cal{L}}_{ij} &=& \sum_{k \neq i,j} f_{C}(r_{ik})~\omega_{ik}~g(\theta_{ijk}) \\
      g(\theta_{ijk}) &=& 1 + \frac{c_{i}^{2}}{d_{i}^{2}} - \frac{c_{i}^{2}}{d_{i}^{2} + (h_{i} - \cos\theta_{ijk})^{2}} \nonumber\end{aligned}

-  **kihs**:

   .. math::

      \begin{aligned}
      \gamma_{ij} &=& (1 + {\cal{L}}_{ij}^{\eta_{i}})^{-\delta_{i}} \nonumber \\
      {\cal{L}}_{ij} &=& \sum_{k \neq i,j} f_{C}(r_{ik})~g(\theta_{ijk})~
      \overbrace{\exp \left[\alpha_{i} (r_{ij}-r_{ik})^{\beta_{i}}\right]}^{\omega_{ik}} \nonumber \\
      g(\theta_{ijk}) &=& c_{1i} + {g}_{o}(\theta_{ijk})~{g}_{a}(\theta_{ijk}) \\
      {g}_{o}(\theta_{ijk}) &=& \frac{c_{2i}~(h_{i} - \cos\theta_{ijk})^{2}}
      {c_{3i} + (h_{i} - \cos\theta_{ijk})^{2}} \nonumber \\
      {g}_{a}(\theta_{ijk}) &=& 1 + {c_{4i}~\exp \left[-c_{5i}~(h_{i} - \cos\theta_{ijk})^{2}\right]}~~, \nonumber\end{aligned}

where the term :math:`{\cal{L}}_{ij}` defines the effective coordination
number of atom :math:`i` i.e. the number of nearest neighbors, taking
into account the relative distance of the two neighbors, :math:`i` and
:math:`k`, :math:`r_{ij}-r_{ik}`, and the bond angle,
:math:`\theta_{ijk}`, between them with respect to the central atom
:math:`i`. The function :math:`g(\theta)` has a minimum for
:math:`h_{i}=\cos(\theta_{ijk})`, the parameter :math:`d_{i}` in
**ters** and :math:`c_{3i}` in **kihs** determines how sharp the
dependence on angle is, whereas the rest express the strength of the
angular effect. Further mixed parameters are defined as:

.. math::

   \begin{aligned}
   a_{ij} = (a_{i} + a_{j})/2&,&~b_{ij} = (b_{i} + b_{j})/2 \nonumber \\
   A_{ij} = (A_{i} A_{j})^{1/2}&,&~B_{ij} = (B_{i} B_{j})^{1/2} \\
   R_{ij} = (R_{i} R_{j})^{1/2}&,&~S_{ij} = (S_{i} S_{j})^{1/2}~~.
   \nonumber\end{aligned}

Singly subscripted parameters, such as :math:`a_{i}` and
:math:`\eta_{i}`, depend only on the type of atom.

For **ters** the chemistry between different atom types is locked in the
two sets bi-atomic parameters :math:`\chi_{ij}` and :math:`\omega_{ij}`:

.. math::

   \begin{aligned}
   \chi_{ii}~=~1&,&~\chi_{ij}~=~\chi_{ji} \nonumber \\
   \omega_{ii}~=~1&,&~\omega_{ij}~=~\omega_{ji}~~,\end{aligned}

which define only one independent parameter each per pair of atom types.
The :math:`\chi` parameter is used to strengthen or weaken the
heteropolar bonds, relative to the value obtained by simple
interpolation. The :math:`\omega` parameter is used to permit greater
flexibility when dealing with more drastically different types of atoms.

The force on an atom :math:`\ell` derived from this potential is
formally calculated with the formula:

.. math::

   f_{\ell}^{\alpha} = -\frac{\partial}{\partial r_{\ell}^{\alpha}}
   E_{\tt tersoff} = \frac{1}{2} \sum_{i}\sum_{j \neq i}
   -\frac{\partial}{\partial r_{\ell}^{\alpha}} U_{ij}~~,

with atomic label :math:`\ell` being one of :math:`i,j,k` and
:math:`\alpha` indicating the :math:`x,y,z` component. The derivative
after the summation is worked out as

.. math::

   -\frac{\partial U_{ij}}{\partial r_{\ell}^{\alpha}} =
   -\frac{\partial}{\partial r_{\ell}^{\alpha}} f_{C}(r_{ij}) f_{R}(r_{ij}) +
    \gamma_{ij} \frac{\partial}{\partial r_{\ell}^{\alpha}} f_{C}(r_{ij}) f_{A}(r_{ij}) +
    f_{C}(r_{ij}) f_{A}(r_{ij}) \frac{\partial}{\partial r_{\ell}^{\alpha}} \gamma_{ij}~~,

with the contributions from the first two terms being:

.. math::

   \begin{aligned}
   -\frac{\partial}{\partial r_{\ell}^{\alpha}} f_{C}(r_{ij}) f_{R}(r_{ij})=&
   -\left\{ f_{C}(r_{ij}) \frac{\partial}{\partial r_{ij}} f_{R}(r_{ij}) +
   f_{R}(r_{ij}) \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) \right\}  \nonumber \\
   &\times ~\left\{ \delta_{j \ell} \frac{r_{i \ell}^{\alpha}}{r_{i \ell}} -
   \delta_{i \ell} \frac{r_{\ell j}^{\alpha}}{r_{\ell j}} \right\}
   \end{aligned}

.. math::

   \begin{aligned}
   \gamma_{ij} \frac{\partial}{\partial r_{\ell}^{\alpha}} f_{C}(r_{ij}) f_{A}(r_{ij})=&
   \gamma_{ij} \left\{ f_{C}(r_{ij}) \frac{\partial}{\partial r_{ij}} f_{A}(r_{ij}) +
   f_{A}(r_{ij}) \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) \right\}  \nonumber \\
   &\times ~\left\{ \delta_{j \ell} \frac{r_{i \ell}^{\alpha}}{r_{i \ell}} -
   \delta_{i \ell} \frac{r_{\ell j}^{\alpha}}{r_{\ell j}} \right\}~~,\end{aligned}

and from the third (angular) term:

-  **ters**:

   .. math::

      \begin{aligned}
      f_{C}(r_{ij}) f_{A}(r_{ij}) \frac{\partial}{\partial r_{\ell}^{\alpha}} \gamma_{ij}=&
      f_{C}(r_{ij}) f_{A}(r_{ij})~\chi_{ij}~~ \nonumber \\
      &\times \left( -\frac{1}{2} \right) \left( 1 + {\beta_{i}}^{\eta_{i}}~{\cal{L}}_{ij}^{\eta_{i}}
      \right)^{-\frac{1}{2 \eta_{i}} - 1} {\beta_{i}}^{\eta_{i}}~{\cal{L}}_{ij}^{\eta_{i}-1}
      \frac{\partial}{\partial r_{\ell}^{\alpha}} {\cal{L}}_{ij}~~,\end{aligned}

   where

   .. math::

      \frac{\partial}{\partial r_{\ell}^{\alpha}} {\cal{L}}_{ij} =
      \frac{\partial}{\partial r_{\ell}^{\alpha}} \sum_{k \neq i,j}
      \omega_{ik}~f_{C}(r_{ik})~g(\theta_{ijk})~~.  \nonumber

   The angular term can have three different contributions depending on
   the index of the particle participating in the interaction:

   .. math::

      \begin{aligned}
      \ell~=~i~&:&~\frac{\partial}{\partial r_{i}^{\alpha}} {\cal{L}}_{ij} = \sum_{k \neq i,j} \omega_{ik}
      \left[ g(\theta_{ijk}) \frac{\partial}{\partial r_{i}^{\alpha}} f_{C}(r_{ik}) +
      f_{C}(r_{ik}) \frac{\partial}{\partial r_{i}^{\alpha}} g(\theta_{ijk}) \right] \nonumber \\
      \ell~=~j~&:&~\frac{\partial}{\partial r_{j}^{\alpha}} {\cal{L}}_{ij} = \sum_{k \neq i,j} \omega_{ik}
      ~f_{C}(r_{ik}) \frac{\partial}{\partial r_{j}^{\alpha}} g(\theta_{ijk}) \\
      \ell~\neq~i,j~&:&~\frac{\partial}{\partial r_{\ell}^{\alpha}} {\cal{L}}_{ij} = \omega_{i \ell}
      \left[ g(\theta_{ij \ell}) \frac{\partial}{\partial r_{\ell}^{\alpha}} f_{C}(r_{i \ell}) +
      f_{C}(r_{i \ell}) \frac{\partial}{\partial r_{\ell}^{\alpha}} g(\theta_{ij \ell}) \right]~~, \nonumber\end{aligned}

-  **kihs**:

   .. math::

      \begin{aligned}
      f_{C}(r_{ij}) f_{A}(r_{ij}) \frac{\partial}{\partial r_{\ell}^{\alpha}} \gamma_{ij}&=&
      f_{C}(r_{ij}) f_{A}(r_{ij})~\times ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \nonumber \\
      & & \left(-\delta_{i}~\eta_{i}\right) \left( 1 + {\cal{L}}_{ij}^{\eta_{i}}\right)^{-\delta_{i} - 1}
      {\cal{L}}_{ij}^{\eta_{i}-1} \frac{\partial}{\partial r_{\ell}^{\alpha}} {\cal{L}}_{ij}~~,\end{aligned}

   where

   .. math::

      \frac{\partial}{\partial r_{\ell}^{\alpha}} {\cal{L}}_{ij} =
      \frac{\partial}{\partial r_{\ell}^{\alpha}} \sum_{k \neq i,j}
      \omega_{ik}(r_{ij},r_{ik})~f_{C}(r_{ik})~g(\theta_{ijk})~~.  \nonumber

   It is worth noting that the derivative of :math:`\omega_{ik}`:

   .. math::

      \frac{\partial}{\partial r_{\ell}^{\alpha}} \omega_{ik} =
      \alpha_{i}~\beta_{i}~(r_{ij}-r_{ik})^{\beta_{i}-1}~\omega_{ik}~
      \left\{ (\delta_{\ell j}-\delta_{\ell i}) \frac{r_{ij}^{\alpha}}{r_{ij}}-
      (\delta_{\ell k}-\delta_{\ell i}) \frac{r_{ik}^{\alpha}}{r_{ik}} \right\}~~,

   now has three different contributions depending on the index of the
   particle participating in the interaction! Hence the angular term’s
   derivative is more elaborate to express than the one in the **ters**
   case.

The derivative of :math:`g(\theta_{ijk})` is worked out in the following
manner:

.. math::

   \frac{\partial}{\partial r_{\ell}^{\alpha}} g(\theta_{ijk}) =
   \frac{\partial g(\theta_{ijk})}{\partial \theta_{ijk}}~
   \frac{-1}{\sin \theta_{ijk}}~\frac{\partial}{\partial r_{\ell}^{\alpha}}
   \left\{ \frac{\underline{r}_{ij} \cdot \underline{r}_{ik}} {r_{ij}~r_{ik}} \right\}~~,

where

.. math::

   \begin{aligned}
   \frac{\partial g(\theta_{ijk})}{\partial \theta_{ijk}}=&
   \frac{2~c_{i}^{2}(h_{i} - \cos \theta_{ijk})~\sin \theta_{ijk}}
   {[d_{i}^{2} + (h_{i} - \cos \theta_{ijk})^{2}]^{2}} \\
   \frac{\partial}{\partial r_{\ell}^{\alpha}}
   \left\{\frac{\underline{r}_{ij}\cdot\underline{r}_{ik}}{r_{ij}r_{ik}}\right\}=&
   (\delta_{\ell j}-\delta_{\ell i})\frac{r_{ik}^{\alpha}}{r_{ij}r_{ik}} +
   (\delta_{\ell k}-\delta_{\ell i})\frac{r_{ij}^{\alpha}}{r_{ij}r_{ik}} \nonumber \\
   & - \cos(\theta_{jik}) \left\{(\delta_{\ell j}-\delta_{\ell i})\frac{r_{ij}^{\alpha}}{r_{ij}^{2}}+
   (\delta_{\ell k}-\delta_{\ell i})\frac{r_{ik}^{\alpha}}{r_{ik}^{2}}\right\}~~.\end{aligned}

The contribution to be added to the atomic virial can be derived as

.. math::

   \begin{aligned}
   {\cal W} = & 3V \frac{\partial E_{\tt tersoff}}{\partial V} =
   \frac{3~V}{2} \sum_{i} \sum_{j \neq i} \frac{\partial U_{ij}}{\partial V} = \sum_{i} \underline{r_{i}} \cdot \underline{f_{i}} \\
   =& \frac{1}{2}\sum_{i} \sum_{j \neq i} -\left( \underline{r_{ij}} \cdot \underline{f_{ij}} + \underline{r_{ik}} \cdot \underline{f_{ik}} \right) = \frac{1}{2}\sum_{i} \sum_{j \neq i} \left(\frac{\partial U_{ij}}{\partial r_{ij}} \cdot \underline{r_{ij}} + \frac{\partial U_{ik}}{\partial r_{ik}} \cdot \underline{r_{ik}}\right) \nonumber\end{aligned}

-  **ters**:

   .. math::

      \begin{aligned}
      {\cal W} =& \frac{1}{2} \sum_{i} \sum_{j \neq i} \left\{ \left[
      \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) f_{R}(r_{ij}) -
      \gamma_{ij} \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) f_{A}(r_{ij}) \right] r_{ij}   \right. \nonumber \\
      & \phantom{xxxxxxxxx} -\left( -\frac{1}{2} \right) f_{C}(r_{ij}) f_{A}(r_{ij})~\chi_{ij}
      \left( 1 + {\beta_{i}}^{\eta_{i}}~{\cal{L}}_{ij}^{\eta_{i}} \right)^{-\frac{1}{2 \eta_{i}} - 1}
      {\beta_{i}}^{\eta_{i}}~{\cal{L}}_{ij}^{\eta_{i}-1}  \\
      & \phantom{xxxxxxxxx}\times \left. \sum_{k \neq i,j} \omega_{ik}~g(\theta_{ijk}) \left[
      \frac{\partial}{\partial r_{ik}} f_{C}(r_{ik}) \right] r_{ik} \right\}~~,
      \nonumber
      \end{aligned}

-  **hiks**:

   .. math::

      \begin{aligned}
      {\cal W} =& \frac{1}{2} \sum_{i} \sum_{j \neq i} \left\{ \left[
      \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) f_{R}(r_{ij}) -
      \gamma_{ij} \frac{\partial}{\partial r_{ij}} f_{C}(r_{ij}) f_{A}(r_{ij}) \right] r_{ij}  \right. \nonumber \\
      & \phantom{xxxxxxxxx} -\left(-\delta_{i}~\eta_{i}\right) f_{C}(r_{ij}) f_{A}(r_{ij})~\chi_{ij}
      \left( 1 + {\cal{L}}_{ij}^{\eta_{i}} \right)^{-\delta_{i} - 1} {\cal{L}}_{ij}^{\eta_{i}-1}  \\
      & \phantom{xxxxxxxxx}\times \left. \sum_{k \neq i,j} \omega_{ik} \left[r_{ik}~g(\theta_{ijk})\frac{\partial}{\partial r_{ik}} f_{C}(r_{ik}) +
      \alpha_{i}~\beta_{i}~(r_{ij}-r_{ik})^{\beta_{i}}f_{C}(r_{ik}) \right] \right\}~~.
      \nonumber
      \end{aligned}

The contribution to be added to the atomic :index:`stress tensor` is given by

.. math:: \sigma^{\alpha \beta} = -r_{i}^{\alpha} f_{i}^{\beta}~~,

where :math:`\alpha` and :math:`\beta` indicate the :math:`x,y,z`
components. The stress tensor is symmetric.

Interpolation arrays, vter and gter (set up in tersoff_generate) -
similar to those in van der Waals interactions (Section `3.1 <#vdw>`__),
are used in the calculation of the Tersoff forces, virial and stress.

.. index:: single: potential;Tersoff

The Tersoff potentials are very short ranged, typically of order
:math:`3` Å. This property, plus the fact that Tersoff potentials (two-
and three-body contributions) scale as :math:`N^{3}`, where :math:`N` is
the number of particles, makes it essential that these terms are
calculated by the link-cell method :cite:`eastwood-80a`.

.. index:: single: potential;Tersoff

DL_POLY_4 applies no long-ranged corrections to the Tersoff potentials.
In DL_POLY_4 Tersoff forces are handled by the routine ``tersoff_forces``.

.. _three-body:

Three-Body Potentials
---------------------

.. index:: single: potential;three-body

The three-body potentials in DL_POLY_4 are mostly valence :index:`angle<potential;valence angle>`
forms. (They are primarily included to permit simulation of amorphous materials
e.g. silicate glasses.) However, these have been extended to include the
:index:`Dreiding<force field;Dreiding>` :cite:`mayo-90a` hydrogen bond. The potential forms
available are as follows:

#. Harmonic: (\ **harm**)

   .. math:: U(\theta_{jik}) = {k \over 2} (\theta_{jik} - \theta_{0})^{2}
      :label: harm_3b_eq

#. Truncated harmonic: (\ **thrm**)

   .. math::
      :label: trunc_harm_3b_eq

      U(\theta_{jik}) = {k \over 2} (\theta_{jik} - \theta_{0})^{2}
      \exp[-(r_{ij}^{8} + r_{ik}^{8}) / \rho^{8}]

#. Screened Harmonic: (\ **shrm**)

   .. math::
      :label: screen_harm_3b_eq

      U(\theta_{jik}) = {k \over 2} (\theta_{jik} - \theta_{0})^{2}
      \exp[-(r_{ij} / \rho_{1} + r_{ik} / \rho_{2})]

#. Screened Vessal :cite:`vessal-94a`: (\ **bvs1**)

   .. math::
      :label: screen_vessal_3b_eq

      \begin{aligned}
      U(\theta_{jik}) =& {k \over 8(\theta_{0}-\pi)^{2}} \left\{ \left[
      (\theta_{0} -\pi)^{2} -(\theta_{jik}-\pi)^{2} \right]^{2} \right\}  \nonumber \\
      & \times\exp[-(r_{ij} / \rho_{1} + r_{ik} / \rho_{2})]\end{aligned}

#. Truncated Vessal :cite:`smith-95a`: (\ **bvs2**)

   .. math::
      :label: trunc_vessal_3b_eq

      \begin{aligned}
      U(\theta_{jik})=& k~\big[ \theta_{jik}^a (\theta_{jik}-\theta_0)^{2}
      (\theta_{jik}+\theta_{0}-2\pi)^{2} \nonumber \\
      & - {a \over 2} \pi^{a-1} (\theta_{jik}-\theta_{0})^{2}(\pi -
      \theta_{0})^{3}\big]~\exp[-(r_{ij}^{8} + r_{ik}^{8}) / \rho^{8}]\end{aligned}

#. :index:`Dreiding<force field;Dreiding>` hydrogen bond :cite:`mayo-90a`: (\ **hbnd**)

   .. math::
      :label: dHBond_3b_eq

      U(\theta_{jik}) =
      D_{hb}~\cos^{4}(\theta_{jik})~[5(R_{hb}/r_{jk})^{12}-6(R_{hb}/r_{jk})^{10}]

.. note:: 
   
   For the hydrogen :index:`bond<potential;chemical bond>`, the hydrogen atom *must* be the
   central atom. 

.. index::
   single: potential;valence angle 
   single: potential;three-body

Several of these functions are identical to those
appearing in the *intra-*\ molecular valence angle descriptions above.
There are significant differences in implementation however, arising
from the fact that the three-body potentials are regarded as
*inter-*\ molecular. Firstly, the atoms involved are defined by atom
types, not specific indices. Secondly, there are *no* excluded atoms
arising from the three-body terms. (The inclusion of other potentials,
for example pair potentials, may in fact be essential to maintain the
structure of the system.)

.. index:: single: potential;three-body

The three-body potentials are very short ranged, typically of order
:math:`3` Å. This property, plus the fact that three-body potentials
scale as :math:`N^{4}`, where :math:`N` is the number of particles,
makes it essential that these terms are calculated by the link-cell
method :cite:`eastwood-80a`.

.. index:: single: potential; valence angle

The calculation of the forces, virial and :index:`stress tensor` as described in
the section valence angle potentials above.

.. index:: single: potential;three-body

DL_POLY_4 applies no long-ranged corrections to the three-body
potentials. The three-body forces are calculated by the routine
``three_body_forces``.

.. _four-body:

Four-Body Potentials
--------------------

.. index:: potential;four-body

The four-body potentials in DL_POLY_4 are entirely inversion 
:index:`angle<potential;inversion>`
forms, primarily included to permit simulation of amorphous materials
(particularly borate glasses). The potential forms available in
DL_POLY_4 are as follows:

#. Harmonic: (\ **harm**)

   .. math:: U(\phi_{ijkn}) = {k \over 2}~(\phi_{ijkn} - \phi_{0})^{2}
      :label: harm_4b_eq

#. Harmonic cosine: (\ **hcos**)

   .. math:: U(\phi_{ijkn}) = {k \over 2}~(\cos(\phi_{ijkn}) - \cos(\phi_{0}))^{2}
      :label: harm_cos_4b_eq

#. Planar potential: (\ **plan**)

   .. math:: U(\phi_{ijkn}) = A~[ 1 - \cos(\phi_{ijkn})]
      :label: planar_4b_eq

These functions are identical to those appearing in the
*intra-*\ :index:`molecular<potential;intramolecular>`
inversion angle descriptions above. There are
significant differences in implementation however, arising from the fact
that the four-body :index:`potentials<potential;four-body>` 
are regarded as *inter-*\ molecular.
Firstly, the atoms involved are defined by atom types, not specific
indices. Secondly, there are *no* excluded atoms arising from the
:index:`four-body<potential;four-body>` terms. (The inclusion of other potentials, for example pair
potentials, may in fact be essential to maintain the structure of the
system.)

.. index:: potential;four-body

The four-body potentials are very short ranged, typically of order
:math:`3~`\ Å. This property, plus the fact that four-body potentials
scale as :math:`N^{4}`, where :math:`N` is the number of particles,
makes it essential that these terms are calculated by the link-cell
method :cite:`eastwood-80a`.

The calculation of the forces, virial and :index:`stress tensor` described in the
section on inversion angle potentials above.

.. index:: single: potential;four-body

DL_POLY_4 applies no long-ranged corrections to the four body
potentials. The four-body forces are calculated by the routine
``four_body_forces``.
