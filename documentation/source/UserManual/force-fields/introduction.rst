Introduction to the DL_POLY_4 Force Field
=========================================

The :index:`force field` is the set of functions needed to define the
interactions in a molecular system. These may have a wide variety of
analytical forms, with some basis in chemical physics, which must be
parameterised to give the correct energy and forces. A huge variety of
forms is possible and for this reason the DL_POLY_4 :index:`force field<force field;DL_POLY>` is
designed to be agnostic and adaptable. While it is not supplied with its
own :index:`force field` parameters, many of the functions familiar to :index:`GROMOS<pair:force field;GROMOS>`
:cite:`gunsteren-87a`, :index:`Dreiding<pair:force field;Dreiding>` :cite:`mayo-90a`
and :index:`AMBER<force field;AMBER>` :cite:`weiner-86a` users have been coded in the
package, as well as less familiar forms. In addition retains the
possibility of the user defining additional potentials.

In DL_POLY_4 the total configuration energy of a molecular system may be
written as:

.. math::
   :label: decomp-ene

   \begin{aligned}
   \label{eq:decomp-ene}
   U(\underline{r}_{1},\underline{r}_{2},\ldots,\underline{r}_{N})=&
         \sum_{i_{shel}=1}^{N_{shel}} U_{shel}(i_{shel},\underline{r}_{core},\underline{r}_{shell}) \nonumber \\
   & + \sum_{i_{teth}=1}^{N_{teth}} U_{teth}(i_{teth},\underline{r}_{i}^{\mathbf{ t}=t},\underline{r}_{i}^{\mathbf{ t}=0}) \nonumber \\
   & + \sum_{i_{bond}=1}^{N_{bond}} U_{bond}(i_{bond},\underline{r}_{a},\underline{r}_{b}) \nonumber \\
   & + \sum_{i_{angl}=1}^{N_{angl}} U_{angl}(i_{angl},\underline{r}_{a},\underline{r}_{b},\underline{r}_{c}) \nonumber \\
   & + \sum_{i_{dihd}=1}^{N_{dihd}} U_{dihd}(i_{dihd},\underline{r}_{a},\underline{r}_{b},\underline{r}_{c},\underline{r}_{d}) \nonumber \\
   & + \sum_{i_{inv}=1}^{N_{inv}} U_{inv}(i_{inv},\underline{r}_{a},\underline{r}_{b},\underline{r}_{c},\underline{r}_{d}) \nonumber \\
   & + \sum_{i=1}^{N-1}\sum_{j>i}^{N} U_{2\textrm{-}body}^{(metal,vdw,electostatics)}(i,j,|\underline{r}_{i}-\underline{r}_{j}|) \\
   & + \sum_{i=1}^{N}\sum_{j{\ne}i}^{N}\sum_{k{\ne}j}^{N} U_{tersoff}(i,j,k,\underline{r}_{i},\underline{r}_{j},\underline{r}_{k}) \nonumber \\
   & + \sum_{i=1}^{N-2}\sum_{j>i}^{N-1}\sum_{k>j}^{N} U_{3\textrm{-}body}(i,j,k,\underline{r}_{i},\underline{r}_{j},\underline{r}_{k}) \nonumber \\
   & + \sum_{i=1}^{N-3}\sum_{j>i}^{N-2}\sum_{k>j}^{N-1}\sum_{n>k}^{N} U_{4\textrm{-}body}(i,j,k,n,\underline{r}_{i},\underline{r}_{j},\underline{r}_{k},\underline{r}_{n}) \nonumber \\
   & + \sum_{i=1}^{N}U_{extn}(i,\underline{r}_{i},\underline{v}_{i})~~,\nonumber\end{aligned}

where
:math:`U_{shel},~U_{teth},~U_{bond},~U_{angl},~U_{dihd},~U_{inv},~U^{(metal)}_{2\textrm{-}body},~U_{tersoff},~U_{3\textrm{-}body}` and :math:`U_{4\textrm{-}body}` are
empirical interaction functions representing ion core-shell
:index:`polarisation<polarisation;shell model>`, tethered :index:`particles<potential;tether>`, chemical :index:`bonds<potential;chemical bond>`, valence :index:`angles<potential;valence angle>`,
:index:`dihedral<potential;dihedral>` (and :index:`improper<potential;improper dihedral>` dihedral angles), inversion :index:`angles<potential;inversion>`, two-body,
:index:`Tersoff<potential;Tersoff>`, :index:`three-body<potential;three-body>` and :index:`four-body<potential;four-body>` forces respectively. The first six are
regarded by DL_POLY_4 as *intra*-molecular interactions and the next
four as *inter*-molecular interactions. The final term :math:`U_{extn}`
represents an *external field* :index:`potential<potential;external field>`. The position vectors
:math:`\underline{r}_{a},\underline{r}_{b},\underline{r}_{c}` and :math:`\underline{r}_{d}`
refer to the positions of the atoms specifically involved in a given
interaction. (Almost universally, it is the *differences* in position
that determine the interaction.) The numbers
:math:`N_{shel},~N_{teth},~N_{bond},~N_{angl}`, :math:`N_{dihd}` and
:math:`N_{inv}` refer to the total numbers of these respective
interactions present in the simulated system, and the indices
:math:`i_{shel},~i_{teth},~i_{bond},~i_{angl},~i_{dihd}` and
:math:`i_{inv}` uniquely specify an individual interaction of each type.
It is important to note that there is no global specification of the
intramolecular interactions in DL_POLY_4- all core-shell units, tethered
particles, chemical bonds, valence angles, dihedral angles and inversion
angles must be individually cited. The same applies for bond :index:`constraints<constraints;bond>`
and PMF :index:`constraints<constraints;PMF>`.


.. index::
      single: potential;non-bonded
      single: potential;van der Waals
      single: potential;electrostatics
      single: potential;metal
      single: potential;three-body
      single: potential;chemical bond

The indices :math:`i`, :math:`j` (and :math:`k`, :math:`n`) appearing in
the intermolecular interactions’ (non-bonded) terms indicate the atoms
involved in the interaction. There is normally a very large number of
these and they are therefore specified globally according to the atom
*types* involved rather than indices. In DL_POLY_4 it is assumed that
the "pure" two-body terms arise from short-ranged interactions such as
van der Waals interactions (or alternatively DPD soft interactions,
coarse-grained interactions, hard-wall nuclear interactions) and
electrostatic interactions (coulombic, also regarded as long-ranged).
Long-ranged forces require special techniques to evaluate accurately
(see Section \ :ref:`coulomb`). The metal terms are many-body
interactions which are functionally presented in an expansion of many
two-body contributions augmented by a function of the local density,
which again is derived from the two-body spatial distribution (and these
are, therefore, evaluated in the two-body routines). In DL_POLY_4 the
three-body terms are restricted to valence angle and H-bond forms.


.. index:: 
      single: force field
      single:  constraints!bond

Throughout this chapter the description of the force field assumes the
simulated system is described as an assembly of atoms. This is for
convenience only, and readers should understand that DL_POLY_4 does
recognize molecular entities, defined through constraint bonds and rigid
bodies. In the case of rigid bodies, the atomic forces are resolved into
molecular forces and torques. These matters are discussed in greater
detail in Sections :ref:`shake-rattle` and
:ref:`rigid`.
