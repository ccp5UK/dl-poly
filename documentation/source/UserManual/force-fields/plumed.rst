.. _plumed:

Free Energy Capabilities via the PLUMED plugin
==============================================

.. index:: single: WWW
   
DL_POLY_4 supports a native integration with PLUMED -
`<http://www.plumed.org/>`_. PLUMED is an open source library for free energy
calculations in molecular systems which works together with some of the
most popular molecular dynamics engines. Free energy calculations can be
performed as a function of many order parameters with a particular focus
on biological problems, using state of the art methods such as
metadynamics, umbrella sampling
:cite:`torrie-77a,kastner-11a` and Jarzynski-equation based
steered MD. The software, written in C++, can be easily interfaced with
both FORTRAN and C/C++ codes.

Using PLUMED can be as simple as adding the keyword **plumed** in your
CONTROL file. By default the input file for PLUMED is called PLUMED and
shall be placed in the same place as your other DL_POLY_4 input files.
Once DL_POLY_4 runs by default OUTPUT.PLUMED will be generated in
addition to the normal PLUMED and DL_POLY_4 output files. The default
names of the files can be changed by using **input** and **log**
parameters with the **plumed** keyword (see
Section :ref:`control-file`).

.. note::
   
   This feature should be considered as **young** rather than
   mature, so please report any bugs and provide feedback regarding its
   improvements.

.. note::
   
   To use the PLUMED functionality within DL_POLY_4 one must
   further ensure that DL_POLY_4 is cross-compiled with it.