.. _equilibration:

Equilibration
=============

Introduction
~~~~~~~~~~~~

DL_POLY_5 includes functionality for equilibration. These include force limiting, temperature (velocity) rescaling and temperature re-gaussing, minimisation (force, energy, or displacement), and progressive heating and cooling. The equilibration statistics may be optionally recorded.

All equilibration functionality applies during the initial equilibration period, with the exception of heating and cooling which may be applied outside of equilibration. This is controlled in the new style CONTROL (see Section :ref:`new-control-file`) by the **time_run** and **time_equilibration** directives. E.g. using **time_equilibration 1000 steps** and **time_run 2000 steps** will simulate in equilibration mode for 1000 steps followed by normal mode for a further 1000 steps.

By default when **time_equilibration** is non-zero only **equilibration_force_cap** will be activated (at 1000 :math:`k_{b}` T / Å). The other options should be manually requested in control.

Force limiting
^^^^^^^^^^^^^^

A global force cap may be applied during equilibration. This option applies a truncation of force magnitudes to the set value. This can be applied using the **equilibration_force_cap** new control directive with units of :math:`k_{b}` T / Å where T is the user temperature unit.

The default force cap is set to 1000.0 :math:`k_{b}` T / Å.

Temperature rescaling, re-gaussing, and zero-k optimisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Temperature rescaling can be applied by setting the **rescale_frequency** to a non-zero number of steps. This will apply a scaling of the atom velocities to match the current instantaneous simulation temperature to the current target temperature during equilibration.

Re-gaussing of the instantaneous system temperature may be performed by random pairwise swaps of the energy scaled momenta of dynamically active particles (excluding massless shells and rigid bodies), by setting **regauss_frequency** to a non-zero number of steps.

zero-k optimisation may be set by a non-zero **reset_temperature_interval** step value. This sets the velocity of a particle (accounting for e.g. rigid bodies and their torques similarly) to:


.. math:: {\bf v} = \begin{cases} {\bf f} \frac{{\bf v}\cdot {\bf f}}{{\bf f}\cdot {\bf f}} & {\bf v}\cdot {\bf f} > 0\\ 0 & \text{otherwise}\end{cases}.
      :label: zero-k

By default all these options are off.

Minimisation
^^^^^^^^^^^^

During equilibration DL_POLY may minimise three properties by the conjugate gradient method. These can be the absolute force, relative energy, or absolute displacement. The criterion can be controlled by **minimistion_criterion** (off, force, energy, or distance) to a given tolerance (in commensurate units) with **minimsation_tolerance**. Minimisation can applied at a certain frequency during equilibration. This can be used as a single step procedure with **time_equilibration 1 steps**.

Minimisation will write out the resulting CONFIG file as CFGMIN after each minimisation pass see Section :ref:`cfgminfile`. This is a CONFIG file containing only position data plus a record of the minimisation process in the header.

By default minimisation is off. The tolerance, step-length, and frequency are 0, -1, and 0 respectively and should be also be set if **minimistion_criterion** is not **off**.

Progressive heating and cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A baseline target temperature is set as **temperature** in Kelvin. A user may include a progressive heating or cooling by using the temperature increment functionality. This allows for a constant temperature increment to be applied beginning at a set timestep (in normal or equilibration mode) and ending at (or below depending on the increment) an increment target temperature. By default the starting step is step 0.

For example given **temperature 11 K** a user can heat from 11 Kelvin in steps of 10 Kelvin every 100 steps, starting at step 100, and ending when reaching 100 K by setting

::

    temperature 11 K
    temperature_increment 10 K
    temperature_increment_frequency 100 steps
    temperature_increment_start 100 steps
    temperature_increment_stop 100 K

This results (given enough timesteps) in a target temperature of 91 K for the rest of the simulation. If the increment was 10 K the final target temperature reached would be 100 K exactly.

The value of the temperature increment should always be positive. For cooling, a higher **temperature** must be supplied with a lower **temperature_increment_stop** value. Negative increments will be converted to positive with a warning issued.

When the temperature increment is applied, the thermostat parameters are updated with the new target temperature immediately. This could be used in conjunction with **rescale_frequency** to also scale velocities on a set schedule.

By default heating/cooling is off and the increment, starting step, and stopping temperature are all 0.