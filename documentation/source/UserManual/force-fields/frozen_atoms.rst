Treatment of Frozen Atoms, Rigid Body and Core-Shell Units
==========================================================

Frozen atoms, core-shell units and rigid body units are treated in a
manner similar to that of the **intra**-molecular interactions due to
their “by site” definition.

DL_POLY_4 allows for atoms to be completely immobilized ( *i.e.*
“frozen” at a fixed point in the MD cell). This is achieved by setting
all forces and velocities associated with that atom to zero during each
MD timestep. Frozen atoms are signalled by assigning an atom a non-zero
value for the freeze parameter in the FIELD file. DL_POLY_4 does not
calculate contributions to the virial or the :index:`stress tensor` arising from
the constraints required to freeze atomic positions. Neither does it
calculate contributions from *intra*- and *inter*- molecular
interactions between frozen atoms. As with the :index:`tethering<potential;tethered>` potential, the
reference position of a frozen site is scaled with the cell vectors in
constant pressure simulations. In the case of frozen rigid bodies, their
“centre of mass” is scaled with the cell vectors in constant pressure
simulations and the positions of their constituent sites are then moved
accordingly.

In DL_POLY_4 the frozen atom option is handled by the subroutine
freeze_atoms.

.. index:: single: equations of motion;Euler
    
The rigid body dynamics (see Section :ref:`rigid`) is resolved
by solving the Eulerian equations of rotational motion. However, their
statics includes calculation of the individual contributions of each
RB’s centre of mass stress and virial due to the action of the resolved
forces on sites/atoms constituting it. These contribute towards the
total system stress and pressure.

.. index:: polarisation;shell model
    
As seen in Section :ref:`shell-models` core-shell units
are dealt with (i) kinetically by the adiabatic shell model or (ii)
statically by the dynamic shell model. Both contribute to the total
system stress (pressure) but in different manner. The former does it via
the kinetic stress (energy) and atomic stress (potential energy) due to
the core-shell spring. The latter via atomic stress (potential energy)
due to the shells move to minimised configuration.