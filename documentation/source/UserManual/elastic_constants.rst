.. _elastic-constants:

Elastic Constants
=================

Introduction
~~~~~~~~~~~~

DL_POLY_5 can utilise its on the fly correlator (see Section :ref:`correlation-functions`) to calculated elastic constants from the stress-fluctuation 
method :cite:`Clavier2017Computation` :cite:`Thompson2022General` :cite:`Pereverzev2022Isothermal`. The
method calculates the elasticity tensor :math:`C_{\alpha\beta\mu\nu}` from the three contributing factors,

.. math:: C_{\alpha\beta\mu\nu} &= \langle C^{B}_{\alpha\beta\mu\nu}\rangle \\&- \frac{V}{k_{\mathrm{B}}T}[\langle \sigma_{\alpha\beta} \sigma_{\alpha\beta}\rangle-\langle \sigma_{\alpha\beta}\rangle\langle \sigma_{\mu\nu}\rangle] \\&+ \frac{2Nk_{\mathrm{B}}T}{V}(\delta_{\alpha\mu} \delta_{\beta\nu}+\delta_{\alpha\nu}\delta_{\beta\mu}).\label{eq-elastic}
   :label: eq-elastic

Where the first term is the Born term, which is defined for pair wise additive potentials :math:`U(r)` as

.. math:: C^{B}_{\alpha\beta\mu\nu} = \frac{1}{V}\sum_{i, j \neq i}\biggl(\frac{\partial^{2} U(r^{ij})}{\partial {r^{ij}}^{2}}-\frac{1}{r^{ij}}\frac{\partial U(r^{ij})}{\partial r^{ij}}\biggr)\frac{r^{ij}_{\alpha} r^{ij}_{\beta} r^{ij}_{\mu} r^{ij}_{\nu}}{{r^{ij}}^{2}}.\label{eq-born}
   :label: eq-born

The second term is the stress-fluctuation term calculated from DL_POLY_5's on the fly correlator. And the final term is the 
Kinetic term. Here :math:`N`, :math:`V` and :math:`T` are system atom count, temperature and volume. :math:`\sigma_{\alpha\beta}` is the 
stress tensor, :math:`r^{ij}_{\alpha}` is the :math:`\alpha` component of inter-atomic separation vector  between atoms :math:`i` and :math:`j`. 
In full generality there may be 21 independent entries in :math:`C_{\alpha\beta\mu\nu}` for an arbitrary system. But for example an isotropic solid with
cubic symmetry will have 3 which can be related to e.g. the shear and bulk modulus.

User Control
~~~~~~~~~~~~

Input
^^^^^

To enable calculation of elastic constants three options are required in the CONTROL (see Section :ref:`new-control-file`).
1. The CONTROL directive **elastic_constants On** must be present. 
2. One or more stress correlation functions must be requested.
3. The CONTROL directive **vdw_method direct** must be present.

Since the Born term, :eq:`eq-born`, requires significant additional computations in DL_POLY_5's Van Der Waals routines, this ensures only absolutely required terms 
are calculated.  

Given these requirements, DL_POLY_5 will automatically calculate the Kinetic and Born terms commensurate with all stress correlations requested,
and reports all the possible elastic constants. E.g. the following CONTROL file snippet will compute the stress-fluctuation term, Born and Kinetic 
terms for only :math:`C_{xxyy}` and :math:`C_{xzxz}`.

::

       correlation_observable [s_xx-s_yy s_xz-xz]
       correlation_block_points [5000 5000]

       elastic_constants On

       vdw_method direct

The computable constants are

.. math:: \begin{matrix}C_{xxxx} & C_{xxyy} & C_{xxzz} & C_{xxyz} & C_{xxzx} & C_{xxxy}\\
                        .        & C_{yyyy} & C_{yyzz} & C_{yyyz} & C_{yyzx} & C_{yyxy}\\
                        .        & .        & C_{zzzz} & C_{zzyz} & C_{zzzx} & C_{zzxy}\\
                        .        & .        & .        & C_{yzyz} & C_{yzzx} & C_{yzxy}\\ 
                        .        & .        & .        & .        & C_{zxzx} & C_{zxxy}\\
                        .        & .        & .        & .        & .        & C_{xyxy}
          \end{matrix}\label{eq-constants}
   :label: eq-constants

Output
^^^^^^

The calculated elastic constants will be written to the COR file. For the example above we may obtain the following 
COR (see Section :ref:`correlation-functions`) snippet.

::
      elasticity_tensor:
            components: [C_xxyy , C_xxzz]
            values: [   10.219979    ,   11.545088    ]
            units: Katm

The elastic constants are computed from equation :eq:`eq-elastic`, and written in, Voigt order.