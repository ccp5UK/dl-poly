
.. _pmf:

.. index:: single: constraints;PMF

Potential of Mean Force (PMF) Constraints and the Evaluation of Free Energy
===========================================================================

A generalization of bond constraints can be made to constrain a system
to some point along a reaction coordinate. A simple example of such a
reaction coordinate would be the distance between two ions in solution.
If a number of simulations are conducted with the system constrained to
different points along the reaction coordinate then the mean constraint
force may be plotted as a function of reaction coordinate and the
function integrated to obtain the free energy for the overall process
:cite:`mccammon-87a`. The PMF constraint force, virial and
contributions to the stress tensor are obtained in a manner analogous to
that for a bond constraint (see previous section). The only difference
is that the constraint is now applied between the centres of two groups
which need not be atoms alone. DL_POLY_4reports the PMF constraint
virial, :math:`{\cal W}_{PMF}`, for each simulation. Users can convert
this to the PMF constraint force from

.. math:: G_{PMF} = \frac {{\cal W}_{PMF}} {d_{PMF}}~~,

where is :math:`d_{PMF}` the constraint distance between the two groups
used to define the reaction coordinate.

The routines ``pmf_shake`` and ``pmf_rattle`` are called to apply
corrections to the atomic positions and respectively the atomic
velocities of all particles constituting PMF units.

.. index::
   single: constraints;bond 
   single: constraints;PMF
   
In presence of both bond constraints and PMF constraints. The constraint
procedures, i.e. SHAKE or RATTLE, for both types of constraints are
applied iteratively in order bonds-PMFs until convergence of
:math:`{\cal W}_{PMF}` reached. The number of iteration cycles is
limited by the same limit as for the bond constraints’ procedures
(SHAKE/RATTLE).