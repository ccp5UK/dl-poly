.. _mdpd:

Many-body Dissipative Particle Dynamics
=======================================

The conservative force in standard DPD depends only upon the species interacting and the interparticle separation, 
which yields a quadratic equation of state. Many-body DPD :cite:`pagonabarraga-01a,trofimov-02a` is a method of providing alternative thermodynamic 
behaviours to DPD particles by making conservative forces additionally depend on local densities. 

The free energy of an inhomogeneous system with density :math:`\rho\left(r\right)` can be defined as the following 
in both the continuous and ensemble-averagd discrete forms:

.. math::
   :label: mdpd_free_energy_eq 

   \mathcal{F} &= \int d\vec{r} \rho(\vec{r}) \psi (\overline{\rho} (\vec{r})) \\
               &= \left\langle \sum_i \psi (\overline{\rho} (\vec{r_i})) \right\rangle

where :math:`\psi\left(\rho\right)` is the free energy per particle in a homogeneous system and 
:math:`\overline{\rho}\left(\vec{r}_{i}\right)` is a function related to the density at (and near) the
position of other particles close to particle :math:`i\left(\tilde{\rho}_{i}\right)` to allow the calculation 
an instantaneous free-energy:

.. math:: 

   \tilde{\mathcal{F}} = \sum_i \psi(\tilde{\rho}_i).

The effective (conservative) force on particle :math:`i` can be obtained from the spatial derivative of the 
free energy, although only the excess part (equivalent to the potential energy :math:`U`) is required since the 
kinetic motion of particles automatically accounts for the ideal contribution: 

.. math:: 

   \vec{F}_i^C = -\frac{\partial \tilde{\mathcal{F}}^{ex} (\left\{ \vec{r}_k \right\})}{\partial \vec{r}_i} = -\sum_j \frac{\partial \psi^{ex} (\tilde{\rho}_j)}{\partial \vec{r}_i}

The force can also be expressed in terms of pairwise interactions, taking a form similar to the standard DPD 
conservativce force :eq:`DPD_eq`:

.. math::
   :label: MDPD_single_eq 

   \vec{F}_{ij}^C = \left( \frac{\partial \psi^{ex} (\tilde{\rho}_i)}{\partial \tilde{\rho}_i} + \frac{\partial \psi^{ex} (\tilde{\rho}_j)}{\partial \tilde{\rho}_j} \right) w^C (r_{ij}) \frac{\vec{r}_{ij}}{r_{ij}} 

The local-density approximation can be defined as a weighted average of instantaneous densities: 

.. math:: 

   \tilde{\rho}_i &= \int d\vec{r} \hspace{0.8mm} w^{\rho} (\mid \vec{r} - \vec{r_i} \mid ) \rho (\vec{r}, \left\{ \vec{r}_k \right\}) \nonumber \\ 
               &= \sum_{j \neq i} \int d\vec{r} \hspace{0.8mm} w^{\rho} (\mid \vec{r} - \vec{r_i} \mid ) \delta (\vec{r} - \vec{r}_j) \nonumber \\ 
   \tilde{\rho}_i &= \sum_{j \neq i} w^{\rho} (r_{ij})

with :math:`w^{\rho}(r)` as the weight function vanishing beyond a cutoff :math:`r_{d}` (which can be equal to or smaller than :math:`r_{c}`) and normalized so that :math:`\int_0^{\infty} 4\pi r^2 w^{\rho} (r) dr = 1`. The 
negative of its derivative with respect to distance :math:`r_{ij}` should also equal the weighting function 
:math:`w^{C}` used in the conservative force given in :eq:`MDPD_single_eq`, i.e. 

.. math:: 

   w^{C} = - \frac{dw^{rho}}{dr_{ij}}.

The most frequently used form for the local density weight function is 

.. math:: 

   w^{\rho} (r_{ij}) = \frac{15r_{c}^{3}}{2 \pi r_d^3} \left( 1 - \frac{r_{ij}}{r_d} \right)^2 ~~~~~ \left(r_{ij} < r_d \right)

which reduces to standard DPD when the excess free energy particle is set to :math:`\psi^{ex} (\tilde{\rho}) = \frac{\pi}{30}A\tilde{\rho}`. This forms the basis of the Two-Parameter model for many-body DPD.

Extending many-body DPD to multiple components is possible by allowing the excess free energies to vary
according to particle species, i.e. :math:`\psi^{ex}_{c}` for component :math:`c`, and generalising 
:eq:`MDPD_single_eq` to include the derivatives for each particle and species in the pair. 

The choice of local densities used to evaluate conservative forces and the component-specific excess free 
energies can be difficult to obtain correctly. The use of component-specific 'partial' local densities, 

.. math:: 
   :label: MDPD_local_density_eq 

   \tilde{\rho}_i^{\alpha} = \sum_{j \in \alpha, j \neq i} w^{\rho} (r_{ij})

can *technically* enable any multiplying parameters for :math:`\psi^{ex}_{c}` and their derivatives to be 
specified for each pair of particle species :cite:`sanyal2018`, i.e. recasting :eq:`MDPD_single_eq` as 

.. math:: 
   :label: MDPD_multiple_eq 

   \vec{F}_{ij}^C = \left( \frac{\partial \psi^{ex}_{c(i)} ( \tilde{\rho}^{\alpha}_{i} )}{\partial \tilde{\rho}_{c(j)}} + \frac{\partial \psi^{ex}_{c(j)} ( \tilde{\rho}^{\alpha}_{j} )}{\partial \tilde{\rho}_{c(i)}} \right) w^C (r_{ij}) \frac{\vec{r}_{ij}}{r_{ij}}

with :math:`c(i)` as the component to which particle :math:`i` belongs (e.g. if :math:`i \in \alpha`, 
:math:`c(i) = \alpha`). However, this approach has been demonstrated to give unphysical results and the use of
different parameters for cross-species interactions can make the forces non-conservative :cite:`warren2013b` (i.e 
the condition :math:`\nabla_{j}\vec{F}_{i} = \nabla_{i}\vec{F}_{j}` would no longer apply).

The use of all-species (total) local densities (:math:`\tilde{\rho}_{i} = \sum_{a}\tilde{\rho}^{a}_{i}`) is 
therefore recommended, giving the following conservative force:

.. math:: 

   \vec{F}_{ij}^C = \left( \frac{\partial \psi^{ex}_{c(i)} ( \left\{ \tilde{\rho}\right\}_i )}{\partial \tilde{\rho}_{j}} + \frac{\partial \psi^{ex}_{c(j)} ( \left\{ \tilde{\rho}\right\}_j )}{\partial \tilde{\rho}_{i}} \right) w^C (r_{ij}) \frac{\vec{r}_{ij}}{r_{ij}}

Specificity among species pairs can still be obtained by varying the weight functions :math:`w^{\rho}` used 
to calculate local densities :cite:`vanya2020`, e.g. by allowing for different :math:`r_{d}` values (:math:`r_{d, ij}`) 
as cutoffs, and the associated conservative weighting functions :math:`w^{C}`. 

Generalised many-body DPD
-------------------------

DL_POLY_5 implements a generalised form of many-body DPD :cite:`vanya2020` where the traditional two-parameter 
many-body DPD model :cite:`warren2003a` is extended to enable the modelling of multiple component systems with 
significant density contrasts, including liquid-solid systems. For the two-parameter model the potential 
(excess free energy) per particle is given by, 

.. math::
   :label: mdpd_two_param_potential_eq 

   \psi^{ex}(\tilde{\rho}) = \frac{\pi r^{4}_{c,ij}}{30}A_{ij}\overline{\overline{\rho}} + 
   \frac{\pi r^{4}_{d}}{30}B\tilde{\rho}^{2}

where :math:`\overline{\overline{\rho}}` is equivalent to :math:`\tilde{\rho}` but with the cutoff set to 
:math:`r_{c,ij}` instead of :math:`r_{d}`. The associated pairwise force is equal to 

.. math:: 
   :label: mdpd_two_param_force_eq 

   \vec{F}^{C}_{ij} = \left[ A_{ij}\left(1 - \frac{r_{ij}}{r_{c,ij}}\right) + B\left(\rho_{i} + \rho_{j}\right)\left(1 - \frac{r_{ij}}{r_{d}}\right)\right]\frac{\vec{r_{ij}}}{r_{ij}}

where :math:`\rho_{i}` and :math:`\rho_{j}` are the local densities considering *all* particle species, 
calculated using a quadratic local density weight function: 

.. math:: 

   w^{\rho}(r_{ij}) = \frac{15}{2\pi r^{3}_{d}}\left(1 - \frac{r_{ij}}{r_{d}}\right)^{2}. 

The terms with :math:`A_{ij}` for both force and potential are calculated as per standard Groot-Warren DPD. By 
setting :math:`A_{ij} < 0` and :math:`B > 0`, a vapour-liquid mixture can be modelled and its equation of state
for a single component is given as 

.. math:: 

   p = \rho k_{B}T + \alpha A\rho^{2} + 2\alpha Br^{4}{d}\left(\rho^{3} - c\rho^{2} + d\right)

where :math:`\alpha \approx 0.101`, :math:`c` and :math:`d` are numerical offsets 
with units of :math:`r_{c}^{-3}` and :math:`r_{c}^{-9}` respectively. 
The value of :math:`d` drops to zero when :math:`A>0`.

The generalised form used in DL_POLY_5 extends this model by; (1) allowing the functional form of the local 
density weighting function to vary with each pair of particle species as a power law 
(:math:`n_{ij} \equiv n_{c(i),c{j}}`) function of distance, 

.. math:: 

   w^{\rho}(r_{ij}) = \frac{(n_{ij} + 1)(n_{ij} + 2)(n_{ij} + 3)r^{3}_{c}}{8\pi r^{3}_{d,ij}} \left(1 - \frac{r_{ij}}{r_{d,ij}}\right)^{n_{ij}},

and (2) choosing power law excess free energies per particle (or wrapping functions) that can vary for each 
individual species (:math:`m_{c(i)}`) i.e. 

.. math:: 

   \psi^{ex}(\tilde{\rho}) = \frac{\pi r^{4}_{c,ij}}{30}A_{ij}\overline{\overline{\rho}} + 
   \frac{B\tilde{\rho}^{m_{c(i)}}}{m_{c(i)}}.

Calcuilating the terms with :math:`A_{ij}` using standard Groot-Warren DPD as before, these lead to the 
following form for the pairwise force:

.. math:: 

   F^{C}_{ij} = \left[ A_{ij} (1 - \frac{r_{ij}}{r_{c,ij}}) + \frac{Bn_{ij}(n_{ij} + 1)(n_{ij} + 2)(n_{ij} + 3)r_{c}^{3}}{8\pi r^{4}_{d,ij}} (\rho^{m_{c(i)} - 1}_{i} + \rho^{m_{c(j)} - 1}_{j}) ( 1 - \frac{r_{ij}}{r_{d,ij}} )^{n_{ij} - 1} \right] \frac{\vec{r_{ij}}}{r_{ij}},


which essentially reduces to the two-parameter form when :math:`n_{ij} = 2` for all pairs of species and 
:math:`m_{c(i)} = m_{c(j)} = 2` for all species. 

This model can be applied to multiple-component systems by specifying :math:`A_{ij}`, :math:`n_{ij}` and 
:math:`r_{d,ij}` for each interacting pair of species and specifying :math:`m_{c(i)}` for each species in 
the ```FIELD``` file. :math:`B` must be set to the same value for all species pairs to ensure the 
many-body forces are conservative, the maximum value specified within the ```FIELD``` file will be used. 

.. note:: 

   The expression for :math:`\psi^{ex}` is this generalised for of many-body DPD differes from the original 
   two-parameter form by a factor of :math:`\frac{\pi r_{d}^{4}}{15r_{c}^{6}}` for the term :math:`B`. Moving to the 
   generalised form from the two-parameter model would require :math:`B` to be multiplied by the same factor. 
   Since DL_POLY_5 only implements the generalised form, care must be taken when using it to apply the two-parameter model, particularly with respect to unit conversions. 

