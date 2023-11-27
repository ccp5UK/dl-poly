.. _ewald-precision:

Choosing Ewald Sum Variables
----------------------------

Ewald sum and SPME
~~~~~~~~~~~~~~~~~~

.. index:: single: Ewald;SPME

This section outlines how to optimise the accuracy of the Smoothed
Particle Mesh Ewald sum parameters for a given simulation..

.. index:: single:: Ewald;optimisation

As a guide to beginners DL_POLY_4 will calculate reasonable parameters
if the **ewald precision** directive is used in the CONTROL file (see
Section :ref:`The CONTROL File<control-file>`). A relative error (see
below) of 10\ :math:`^{-6}` is normally sufficient so the directive

:: 
   
   ewald precision 1d-6

will make DL_POLY_4 evaluate its best guess at the Ewald parameters
:math:`\alpha`, ``kmaxa``, ``kmaxb`` and ``kmaxc``, or their doubles if
**ewald** rather than **spme** is specified. (The user should note that
this represents an *estimate*, and there are sometimes circumstances
where the estimate can be improved upon. This is especially the case
when the system contains a strong directional anisotropy, such as a
surface.) These four parameters may also be set explicitly by the
**ewald sum** directive in the CONTROL file. For example the directive

:: 
   
   ewald sum 0.35 6 6 8

which is equvalent to

:: 
   
   spme sum 0.35 12 12 16

would set :math:`\alpha=0.35` Å\ :math:`^{-1}`, :math:`\texttt{kmaxa}=12`,
:math:`{\tt
kmaxb}=12` and :math:`\texttt{kmaxc}=16`\  [1]_. The quickest check on
the accuracy of the Ewald :index:`sum<Ewald;summation>` is to compare the coulombic energy
(:math:`U`) and virial (:math:`\cal W`) in a short simulation. Adherence
to the relationship :math:`U=-{\cal W}`, shows the extent to which the
Ewald :index:`sum<Ewald;summation>` is correctly converged. These variables can be found under the
columns headed ``eng_cou`` and ``vir_cou`` in the OUTPUT file (see
Section :ref:`The OUTPUT FILES<output-files>`).

The remainder of this section explains the meanings of these parameters
and how they can be chosen. The Ewald :index:`sum<Ewald;summation>` can only be used in a three
dimensional periodic system. There are five variables that control the
accuracy: :math:`\alpha`, the Ewald convergence parameter;
:math:`r_{\rm cut}` the real space force cutoff; and the ``kmaxa``,
``kmaxb`` and ``kmaxc`` integers that specify the dimensions of the SPME
charge array (as well as FFT arrays). The three integers effectively
define the range of the reciprocal space sum (one integer for each of
the three axis directions). These variables are not independent, and it
is usual to regard one of them as pre-determined and adjust the others
accordingly. In this treatment we assume that :math:`r_{\rm cut}`
(defined by the **cutoff** directive in the CONTROL file) is fixed for
the given system.

The Ewald sum splits the (electrostatic) sum for the infinite, periodic,
system into a damped real space sum and a reciprocal space sum. The rate
of convergence of both sums is governed by :math:`\alpha`. Evaluation of
the real space sum is truncated at :math:`r=r_{\rm cut}` so it is
important that :math:`\alpha` be chosen so that contributions to the
real space sum are negligible for terms with :math:`r>r_{\rm cut}`. The
relative error (:math:`\epsilon`) in the real space sum truncated at
:math:`r_{\rm cut}` is given approximately :index:`by<Ewald;optimisation>`

.. math::
   :label: relerr_eq

   \epsilon \approx {\rm erfc}(\alpha~r_{\rm cut})/r_{\rm cut}
   \approx \exp[-(\alpha~r_{\rm cut})^{2}]/r_{\rm cut}~~, \label{relerr}

which reciprocally gives an estimate for :math:`\alpha` for a given
:math:`\epsilon`:

.. math:: \alpha \approx \frac{\sqrt{|{\rm ln}(\epsilon~r_{\rm cut})|}}{r_{\rm cut}}~~.

The recommended value for :math:`\alpha` is :math:`3.2/r_{\rm cut}` or
greater (too large a value will make the reciprocal space sum very
slowly convergent). This gives a relative error in the energy of no
greater than :math:`\epsilon = 4 \times 10^{-5}` in the real space sum.
When using the directive **ewald precision** DL_POLY_4 makes use of a
more sophisticated approximation:

.. math:: {\rm erfc}(x) \approx 0.56 \; \exp(-x^{2})/x

to solve recursively for :math:`\alpha`, using
equation :eq:`relerr_eq` to give the first guess.

The relative error in the reciprocal space term is approximately

.. math:: \epsilon \approx \exp(- k_{max}^{2}/4\alpha^{2})/k_{max}^{2}

where

.. math:: k_{max} = \frac{2\pi}{L}~\frac\texttt{kmax}{2}

is largest :math:`k`-vector considered in reciprocal space, :math:`L` is
the width of the cell in the specified direction and ``kmax`` is an
integer.

For a relative error of :math:`4 \times 10^{-5}` this means using
:math:`k_{max} \approx 6.2~\alpha`. ``kmax`` is then

.. math:: \texttt{kmax} > 6.4~L/r_{\rm cut}.

In a cubic system, :math:`r_{\rm cut}=L/2` implies
:math:`\texttt{kmax}=14`. In practice the above equation slightly over
estimates the value of ``kmax`` required, so optimal values need to be
found experimentally. In the above example :math:`\texttt{kmax}=10` or
:math:`12` would be adequate.

If you wish to set the Ewald parameters manually (via the **ewald sum**
or **spme sum** directives) the recommended approach is as follows.
Preselect the value of :math:`r_{\rm cut}`, choose a working a value of
:math:`\alpha` of about 3.2/\ :math:`r_{\rm cut}` and a large value for
the ``kmax`` (say 20 20 20 or more). Then do a series of ten or so
*single* step simulations with your initial configuration and with
:math:`\alpha` ranging over the value you have chosen plus and minus
20%. Plot the Coulombic energy (-:math:`\cal W`) versus :math:`\alpha`.
If the Ewald sum is correctly converged you will see a plateau in the
plot. Divergence from the plateau at small :math:`\alpha` is due to
non-convergence in the real space sum. Divergence from the plateau at
large :math:`\alpha` is due to non-convergence of the reciprocal space
sum. Redo the series of calculations using smaller ``kmax`` values. The
optimum values for ``kmax`` are the smallest values that reproduce the
correct Coulombic energy (the plateau value) and virial at the value of
:math:`\alpha` to be used in the simulation. Note that one needs to
specify the three integers (``kmaxa``, ``kmaxb``, ``kmaxc``) referring
to the three spatial directions, to ensure the reciprocal space sum is
equally accurate in all directions. The values of ``kmaxa``, ``kmaxb``
and ``kmaxc`` must be commensurate with the cell geometry to ensure the
same minimum wavelength is used in all directions. For a cubic cell set
``kmaxa`` = ``kmaxb`` = ``kmaxc``. However, for example, in a cell with
dimensions 2A = 2B = C, (ie. a tetragonal cell, longer in the c
direction than the a and b directions) use 2\ ``kmaxa`` = 2\ ``kmaxb`` =
``kmaxc``.

If the values for the kmax used are too small, the Ewald sum will
produce spurious results. If values that are too large are used, the
results will be correct but the calculation will consume unnecessary
amounts of cpu time. The amount of cpu time increases proportionally to
:math:`\texttt{kmaxa}  \times \texttt{kmaxb} \times \texttt{kmaxc}`.

It is worth noting that the working values of the k-vectors may be
larger than their original values depending on the actual processor
decomposition. This is to satisfy the requirement that the k-vector/FFT
transform down each direction per domain is a multiple of 2, 3 and 5
only, which is due to the GPFA code (single 1D FFT) which the DaFT
implementation relies on. This allowes for greater flexiblity than the
power of 2 multiple restriction in DL_POLY_4 predicessor, DL_POLY_3. As
a consequence, however, execution on different processor decompositions
may lead to different working lengths of the k-vectors/FFT transforms
and therefore slightly different SPME forces/energies whithin the same
level of SPME/Ewald precision/accuracy specified. 

.. note:: 
   
   Although the number of processors along a dimension of the DD grid may be any
   number, numbers that have a large prime as a factor will lead to
   inefficient performance!

.. [1]
   **Important note**: As the SPME method substitues the standard Ewald
   the values of ``kmaxa``, ``kmaxb`` and ``kmaxc`` are the double of
   those in the prescription of the standard Ewald since they specify
   the sides of a cube, not a radius of convergence.
