.. _coarse-graining:

Course Graining
===============

.. _cg-intro:

User-Defined Coarse-Grain Models with Tabulated Force-Fields
------------------------------------------------------------

One can use DL_POLY_4 for preparing and running simulations of
(numerically) coarse-grained (CG) models by using tabulated effective
force-fields (FF) derived from either **(i)** potentials of mean force
(PMF); or **(ii)** iteratively optimised CG models.

In outline, systematic coarse-graining (SCG) of an atomistic system
implies the application of a geometrical projection, or “mapping”, of
the original system onto a considerably reduced set of degrees of
freedom (DoF), whence referred to as a “coarse-grained” model. The
procedure is to recover the average configurational (often termed
“physical”) and topological (often termed “chemical”) force-field
related properties of the original full-atom (FA) model. Integrating out
degrees of freedom ultimately leads to loss of information. However, if
this is done selectively and consistently, the intrinsic information
pertaining to the dominating interactions within the FA system (that
govern its behaviour towards phenomena of interest) and the most
important thermodynamic properties are retained. Then the greatly
reduced phase-space of the CG model system allows for efficient access
to much longer length and time scales for investigating microscopic
phenomena using of the well-established classical MD machinery. The
reduced DoF mapping leads to the generation of an effective CG FF in
terms of the effective interaction potentials and forces between the CG
particles, which are often tabulated numerically.

.. index:: single:: WWW

The initial coarse-grain mapping of the original FA trajectory can be
done with the aid of DL_CGMAP tool – `<http://www.ccp5.ac.uk/software/>`_.
After the CG mapped trajectory has been obtained, the relevant
distribution analysis and the Boltzmann Inversion procedure (producing
tabulated PMF:s) can be performed by using in “replay history” mode,
with the corresponding directives in CONTROL file. For further iterative
optimisation of the CG model the user is advised to use DL_POLY_4 as
simulation engine within the framework of VOTCA package –
`<http://www.votca.org/>`_, which provides a handful of various systematic
coarse-graining methodologies. For more details the user should refer to
the manuals of these tools.

This section describes how to use DL_POLY_4 for two SCG tasks:

-  Post-simulation analysis of the intramolecular (bonded and angular)
   mean-force interactions, based on the calculation of the
   corresponding intramolecular distributions.

-  Preparation and use of the tabulated intramolecular potentials.

.. note::
   
   Although these steps are also applicable to atomistic systems,
   which can be useful for benchmarking purposes, below we shall assume
   that the CG mapping has been done and the following data files have been
   generated for the CG mapped model system: CONTROL, FIELD, CONFIG, and,
   alternatively, HISTORY.

.. _IPDF-analysis:

Intramolecular Probability Distribution Function (PDF) Analysis
---------------------------------------------------------------

Albeit the distribution analysis can be performed on the fly, in the
course of a simulation run, for CG purposes it has to be done on an
existing CG mapped trajectory, by invoking the **replay history** and
**analysis** directives, see
Section :ref:`control-file`. To trigger the PDF
collection and subsequent output of PMF:s for any of the following
intramolecular interactions: bonds, angles, dihedrals and/or inversions,
the CONTROL file must contain any combination of the options
demonstrated in the example below. 

.. note::
   
   Such analyses will only
   be carried out if the desired intramolecular types of interactions are
   defined within FIELD!

::

   TITLE: EXAMPLE OF DL_POLY_4  PDF ANALYSIS DIRECTIVES SNIPPET

   # DIRECTIVES TO INVOKE INTRAMOLECULAR PDF ANALYSIS BY TYPE
   analyse  bonds       sample every 100  nbins 250  rmax 5.0
   analyse  angles      sample every 100  nbins 360  # [ 0 : pi]
   analyse  dihedrals   sample every 100  nbins 720  # [-pi: pi]
   analyse  inversions  sample every 100  nbins 360  # [ 0 : pi]

   # DIRECTIVES TO INVOKE INTRAMOLECULAR PDF ANALYSIS FOR ALL TYPES
   analyse  all         sample every 100  nbins 1000  rmax 5.0

   # DIRECTIVES TO INVOKE PRINTING FOR ANY DEFINED
   # INTER-(RDF->VDW) & INTRA-(bonded) MOLECULAR PDF ANALYSIS
   print analysis

The **analyse** directive acts in a similar manner to that of the
**rdf** directive outlining **(i)** the frequency of sampling (in
steps/frames) and **(ii)** the number of bins over the cutoff interval
of the specified interaction. It is only for bonded pairs that this
cutoff needs specifying, whereas for angles the possible ranges are
known *a priory* by definition. If no cutoff is supplied for bonds, then
it defaults to 2 Å. It is worth noting that, for the sake of accuracy,
the number of bins per PDF must be larger than or equal to
:math:`N_{\rm min} = {\tt Nint} \left( \frac{\rm PDF~cutoff}{\Delta_{\rm max}} \right)`,
where the max bin size, :math:`\Delta_{\rm max}`, is defined internally
for each type of PDF (see setup_module.f90). Otherwise, it will default
to :math:`{\tt max}(N_{\rm min},N_{\rm tab})`, where :math:`N_{\rm tab}`
is the number of bins on the grid defined in the corresponding tabulated
force-field data file (TAB*), if it is provided, otherwise
:math:`N_{\rm tab}=0`. In particular, for bonds :math:`\Delta_{\rm max}`
is 0.01 Å and for angles it is :math:`0.2 \pi/180`, i.e. 0.2 degree.

If **analyse all** option is used in conjunction with any specific
directive, then it triggers the analysis on all PDFs while enforcing on
the individually targeted PDF(s) the following parameters: **(i)** its
sampling interval (frames) – only if it is smaller than, and **(ii)**
its grid number – only if it is larger than those specified for the
targeted PDFs. In the case of bonds it will also enforce the grid range
**rmax** – only if it is larger than that specified in the **analyse
bonds** directive. Hence, the **analyse all** directive allows to
quickly override and/or unify the sampling frequencies and max grid
numbers for all PDFs, provided its parameters improve on the accuracy of
the collected data (compared to the specifications for the individually
targeted PDF:s).

While the statistics are always collected and stored for the future use
in the (binary) REVIVE file for all targeted PDF:s (see
Section :ref:`revive-file`), using the **print
analysis** directive will instruct DL_POLY_4 to additionally print the
data in the OUTPUT and \*DAT files, the latter containing each type of
PDF:s separately (the asterisk stands for one of the following: BND,
ANG, DIH, or INV). As a result, apart from OUTPUT, three data files will
be created for each of the targeted distribution types: BNDDAT,
BNDPMF & BNDTAB – for bonds, ANGDAT, ANGPMF & ANGTAB – for angles,
DIHDAT, DIHPMF & DIHTAB – for dihedrals, and INVDAT, INVPMF & INVTAB –
for inversions; the \*PMF and \*TAB files containing tabulated data for
the respective potential of mean force and force/virial (*TAB files
bearing the data from \*PMF files but *resampled onto a finer grid*; see
the last paragraph in this section).

Partial examples of the \*DAT files for bonds and angles are given
below.

::

   [user@host]$ more BNDDAT
   # TITLE: Hexane FA OPLSAA -> CG mapped with 3 beads (A-B-A)
   # BONDS: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dr_bin
   # bins, cutoff, frames, types:        250    5.000       2285          1
   #
   # r(Angstroms)  PDF_norm(r)  PDF_norm(r)/dVol(r)   @   dr_bin =  0.02000
   #

   # type, index, instances: A        B                 1       2000
       0.01000  0.000000E+00  0.000000E+00
       0.03000  0.000000E+00  0.000000E+00
       0.05000  0.000000E+00  0.000000E+00
   ...
       4.95000  0.000000E+00  0.000000E+00
       4.97000  0.000000E+00  0.000000E+00
       4.99000  0.000000E+00  0.000000E+00

   [user@host]$ more ANGDAT
   # TITLE: Hexane FA OPLSAA -> CG mapped with 3 beads (A-B-A)
   # ANGLES: Probability Density Functions (PDF) := histogram(bin)/hist_sum(bins)/dTheta_bin
   # bins, cutoff, frames, types:        360        180       2285          1
   #
   # Theta(degrees)  PDF_norm(Theta)  PDF_norm(Theta)/sin(Theta)   @   dTheta_bin =  0.50000
   #

   # type, index, instances: A        B        A                 1       1000
       0.25000  0.000000E+00  0.000000E+00
       0.75000  0.000000E+00  0.000000E+00
       1.25000  0.000000E+00  0.000000E+00
   ...
     178.75000  1.380569E-02  6.328564E-01
     179.25000  8.368490E-03  6.393238E-01
     179.75000  2.901532E-03  6.649842E-01

One can see that all the header lines are commented out due to starting
with the hash symbol, **#**. Nonetheless, the header contains some
useful information. The title is, as usual, placed in the first line,
which is followed by an explanatory line with the definition of a
normalised PDF. The third line provides the four most important
descriptors: the number of *bins* on the histogram grid, the *cutoff*
interval (absolute value of the span) over which the distributions are
sampled, the number of *frames* (samples) used, and the number of unique
unit *types* analysed (where “unit” is one of the following: bonds,
angles, dihedrals or inversions). The last explanatory line in the
header, found in between two empty commented-out lines, defines the
meaning of the columns and, at the end, the grid bin size.

The PDF histograms within a \*DAT file are separated by uncommented
empty lines, which makes it possible to directly import and plot all the
data as separate lines in [Xm]Grace 2D plotter. For clarity, the data of
each histogram are preceded by a commented-out line specifying the
*type*, *index* and number of *instances* of the analysed interaction
unit (the latter being counted over the entire system).

In all \*DAT files the first column bears the *bin-centered* abscissa
values (distance or angle), the second column is the distribution
histogram normalized to unity, i.e. its *integral* equals 1, whereas the
third column, if any, contains the PDF data corrected for the volumetric
(entropic) degeneracy of the grid points. Thus, it is the data of the
last column found that are used for calculating PMF
:math:`\sim -\ln`\ (PDF).

The OUTPUT file will contain copies of the first and second columns of
all collected PDF:s, but normalised so that the figures in the second
column *sum up* to 1, which can be checked by examining the third column
as it bears the running sum of a PDF histogram.

In addition to the PDF:s, DL_POLY_4 also calculates the respective PMF:s
and pairwise force dependencies (virial for bonds), which are stored in
the \*PMF and, upon resampling onto a *bin-edge* grid, \*TAB files. **It
is important to note** that, unlike the \*PMF files containing the bare
:math:`-\ln`\ (PDF) data (converted to the requested energy units) *on
the same grid as the PDF histograms*, the force-field data in \*TAB
files are *resampled* onto :math:`\max(N_{\rm min},N_{\rm tab})` grid
points located at the bin edges (as expected by DL_POLY_4 when reading
the potential and force tables), where :math:`N_{\rm tab}` is the grid
number for the respective intramolecular unit type read-in from its
TAB\* file, if provided (otherwise :math:`N_{\rm tab}=0`). Thus, the
\*TAB files obey the DL_POLY_4 format for numerically defined
intramolecular force-field tables (TAB*, see below) and, hence, can be
directly used as such upon renaming: BNDTAB \ :math:`\to` TABBND,
ANGTAB \ :math:`\to` TABANG, DIHTAB \ :math:`\to` TABDIH and
INVTAB \ :math:`\to` TABINV. **The user is, however, strongly advised to
check the quality of the obtained tabulated force-fields before using
those as input for a CG simulation.** Albeit DL_POLY_4 uses a simple
smoothing algorithm for noise-reduction in PMF:s and implements capping
of the forces in the regions of zero-valued PDF:s, in undersampled
regions the PDF and PMF data are likely to suffer from inaccuracy and
increased noise which, most often, require extra attention and
re-fitting manually.

The general format of the above discussed files is shown in
Section :ref:`bonded-files`.

.. _bonded-tables:

Setting up Tabulated Intramolecular Force-Field Files
-----------------------------------------------------

For a user-defined, e.g. coarse-grained, model system the effective
potentials must be provided in a tabulated form. For non-bonded
short-range (VdW) interactions the TABLE file must be prepared as
described in Section :ref:`table-file`. However, the
tabulated data format for intramolecular interactions (bonds, bending
angles, dihedral and inversion angles in a polymer) differs from that of
the TABLE file and assumes three columns: abscissa (distance in Å or
angle in degrees), and two ordinates: potential and force data (virial
for distance dependent interactions – e.g. bonds, and force for angle
dependent interactions – e.g. angles). Shown below are examples of
TABBND and TABANG files, corresponding to the above PDF examples.


.. note::
   
   The PMF and force data have been resampled onto a finer
   grid with points located *at bin edges*.

::

   [user@host]$ more TABBND
   # TITLE: Hexane FA OPLSAA -> CG mapped with 3 beads (A-B-A)
   # 5.0 500

   # A B
    1.00000e-02  -9.0906600e+02   1.3954000e+00
    2.00000e-02  -9.1046200e+02   2.7908000e+00
    3.00000e-02  -9.1185700e+02   4.1862000e+00
   ...

   [user@host]$ more TABANG
   # TITLE: Hexane FA OPLSAA -> CG mapped with 3 beads (A-B-A)
   # 1000

   # A B A
    1.80000e-01   8.8720627e+01   6.9119576e-01
    3.60000e-01   8.8596227e+01   6.9119227e-01
    5.40000e-01   8.8471827e+01   6.9118704e-01
   ...

The input tables for bonds, angles, dihedrals and inversions are named
TABBND, TABANG, TABDIH and TABINV, correspondingly. The format of these
files is fixed in terms of the line, or *record*, order. In particular,
the initial two header lines must contain a title and a record with the
grid specification, and each of the following blocks of tabulated data
must be preceded by an empty line and a one-line descriptor record
containing the (white-space delimited) names of the atoms making up the
given intermolecular interaction unit (same as a *unit type* in \*DAT
files). These descriptor lines can be commented out or not, i.e. having
**#** as the first symbol would not affect the reading operation, but it
would ease importing of the data for plotting and manipulating in
[Xm]Grace software.

In all TAB\* files the number of grid points (bins) must be specified on
the second line (commented-out or not). For angles (TABANG, TABDIH,
TABINV) no other information needs to be provided as their ranges are
pre-determined: :math:`0 < {\Theta\ \rm in\ TABANG\ \&\ TABINV} \le 180`
and :math:`-180 < {\Theta\ \rm in\ TABDIH} \le 180`. In the TABBND file,
however, the “bond cutoff”, :math:`r_{\rm max}` (Å), must be precede the
grid number.

.. note::
   
   All potential and force data are to be provided in the
   same energy units as specified by the user in the FIELD file (see
   Section :ref:`field-file`) with distances in Å and angles
   in degrees. All the data related to angles are internally transformed
   and handled by DL_POLY_4 with angles measured in *radians*.

Finally, in order to instruct DL_POLY_4 to use tabulated intramolecular
force-fields read from the TAB\* files the user has to specify in the
FIELD file the keyword **tab** or **-tab** for each intramolecular
interaction thereby chosen for tabulation (similarly to how it is
described in Section :ref:`field-file`). The dash symbol
(-) in front of the keyword **tab** is only valid for bonds and angles,
and is interpreted in the same manner as in
Table :numref:`(%s)<bond-table>` and
Table :numref:`(%s)<angle-table>`.

.. note::
   
   The VOTCA package is also capable of collecting both intra-
   and inter-molecular stats and producing correct TAB\* files, provided
   the FIELD and HISTORY files exist (albeit VOTCA saves the distributions
   in a format different from \*DAT files).

Below we summarise the sequence of operations the user has to follow in
order to perform the CG distribution analysis and prepare the TAB\*
files for a newly coarse-grained system.

-  Perform CG mapping of the original FA system with the aid of DL_CGMAP
   or VOTCA (in the case of using VOTCA follow its manual; the remainder
   of the list describes using DL_POLY_4 only);

-  Move the data files for the newly CG-mapped system (FIELD_CG,
   CONFIG_CG, HISTORY_CG) into a separate directory under the standard
   names (FIELD, CONFIG, HISTORY) and create the corresponding CONTROL
   file containing the **analysis** and **replay history** directives.

-  As no TAB\* files are yet available for the CG system at this stage,
   one can either **(i)** create the initial TAB\* files padded with
   zeros, or **(ii)** use a FIELD file with fictitious records for all
   the interactions to be tabulated, with the interaction keywords and
   parameters chosen arbitrarily in accord with
   Table :numref:`(%s)<bond-table>`,
   Table :numref:`(%s)<angle-table>`,
   Table :numref:`(%s)<dihedral-table>` and
   Table :numref:`(%s)<inversion-table>`. 

   .. note::

      DL_CGMAP (as well as VOTCA) creates FIELD_CG files that already
      contain interaction descriptors with the **tab** keyword(s) in place,
      so if the route *(ii)* is chosen, the user needs to replace those
      records with the fictitious ones (it is advisory to store the initial
      FIELD_CG for the future use, when the actual TAB\* files are ready).

-  Run DL_POLY_4 with the **replay history**, **rdf** and/or **analysis**
   options invoked in the CONTROL file, which will result in creation of
   the targeted inter- and intra-molecular PDF data files (RDFDAT,
   BNDDAT, ANGDAT, DIHDAT, INVDAT) and the respective PMF files as
   described above. **Note** that only and only when the **rdf** and
   **analysis** options are both active then the VDWPMF and VDWTAB files
   (derived from RDF:s) will be produced, along with RDFDAT. They are
   structured in the same manner and format as their intramolecular
   counterparts. The user can then convert the VDWTAB file into a
   correctly formatted TABLE file by using the utility called pmf2tab.f
   (subject to compilation; found in DL_POLY_4 directory utility) as
   follows.

   ::

      [user@host]$ pmf2tab.exe < VDWTAB

-  Check the data for accuracy and amend the tabulated force-fields.
   Redo the analysis on coarser/finer grid(s), if necessary.

-  When satisfied with the created TAB\* files, run a
   DL_POLY_4 simulation for the prepared CG (or simply user-defined)
   model system.
