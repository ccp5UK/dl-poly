.. _kim:

Open Knowledgebase of Interatomic Models - OpenKIM
==================================================

.. index:: WWW

DL_POLY_4 force-field allows for model interactions specification by
using an OpenKIM model -
`<https://openkim.org/>`_. A KIM model contains all necessary (non-bonded)
interactions and their parameters for a specific model within a
designated container. Thus a KIM model can be used as a force-field
container when made available to DL_POLY_4 at run time (see the
description of the FIELD file in
Section :ref:`field-file`) upon a specification within
FIELD together with a matching molecular description for the model
system, also specified in FIELD. Due to the history of the OpenKIM
initiative and the constraints of the logic of the DL_POLY_4 FIELD
file the designated place for a KIM model specification is in the
non-bonded interactions section :ref:`non-bonded_sec`.

Employing OpenKIM interatomic models (IMs) provides DL_POLY_4 users with
multiple benefits, including:

Reliability
-----------

-  All content archived in OpenKIM is reviewed by the KIM Editor for
   quality.

-  IMs in OpenKIM are archived with full provenance control. Each is
   associated with a maintainer responsible for the integrity of the
   content. All changes are tracked and recorded.

-  IMs in OpenKIM are exhaustively tested using KIM Tests that compute a
   host of material properties, and KIM Verification Checks that provide
   the user with information on various aspects of the IM behavior and
   coding correctness. This information is displayed on the IM’s page
   accessible through the OpenKIM browse interface.

Reproducibility
---------------

-  Each IM in OpenKIM is issued a unique identifier (KIM ID), which
   includes a version number (last three digits). Any changes that can
   result in different numerical values lead to a version increment in
   the KIM ID. This makes it possible to reproduce simulations since the
   specific version of a specific IM used can be retrieved using its KIM
   ID.

-  OpenKIM is a member organization of DataCite and issues digital
   object identifiers (DOIs) to all IMs archived in OpenKIM. This makes
   it possible to cite the IM code used in a simulation in a
   publications to give credit to the developers and further facilitate
   reproducibility.

Currently, DL_POLY_4 supports one type of IM archived in OpenKIM, which
is called a KIM Portable Model (PM). A KIM PM is an independent computer
implementation of an IM written in one of the languages supported by KIM
(C, C++, Fortran) that conforms to the KIM Application Programming
Interface (KIM API) Portable Model Interface (PMI) standard. A KIM PM
will work seamlessly with any simulation code that supports the KIM
API/PMI standard (including ; :index:`see<WWW>` `complete list of supported
codes <https://openkim.org/projects-using-kim/>`__).

OpenKIM IMs are uniquely identified by a KIM ID. The extended KIM ID
consists of a human-readable prefix identifying the type of IM,
authors, publication year, and supported species, separated by two
underscores from the KIM ID itself, which begins with an IM code (MO
for a KIM Portable Model) followed by a unique 12-digit code and a
3-digit version identifier
(e.g. **SNAP_ChenDengTran_2017_Mo__MO_698578166685_000**).

Each OpenKIM IM has a dedicated “Model Page” on OpenKIM providing all
the information on the IM including a title, description, authorship
and citation information, test and verification check results,
visualizations of results, a wiki with documentation and user
comments, and access to raw files, and other information. The URL for
the Model Page is constructed from the extended KIM ID of the IM
:index:`as<WWW>` `<https://openkim.org/id/extended_KIM_ID>`_.

For example, for the spectral neighbor analysis potential (SNAP)
listed above the Model Page is located :index:`at<WWW>`:
`<https://openkim.org/id/SNAP_ChenDengTran_2017_Mo__MO_698578166685_000>`_.

See the `current list of KIM PMs archived in
OpenKIM <https://openkim.org/browse/models/by-species>`_. This list is
sorted by species and can be filtered to display only IMs for certain
species combinations. You can also see `Obtaining KIM
Models <https://openkim.org/doc/usage/obtaining-models/>`_ to learn how
to install a pre-build binary of the OpenKIM repository of models.

.. note::
   It is also possible to locally install IMs not archived in
   OpenKIM, in which case their names do not have to conform to the KIM ID
   format.

.. index:: single: WWW

To download, build, and install the KIM library on your system, see the
`detailed
instructions <https://openkim.org/doc/usage/obtaining-models/>`_ .

.. note::
   
   To use OpenKIM functionality within DL_POLY_4 one must
   further ensure that DL_POLY_4 is compiled with OpenKIM support (see
   building.md).

Citation of OpenKIM IMs
-----------------------

.. index:: single: WWW
   
When publishing results obtained using OpenKIM IMs researchers are
requested to cite the OpenKIM project
`(Tadmor) <https://link.springer.com/article/10.1007%2Fs11837-011-0102-6>`__,
KIM API `(Elliott) <https://doi.org/10.25950/FF8F563A>`__, and the
specific IM codes used in the simulations, in addition to the relevant
scientific references for the IM. The citation format for an IM is
displayed on its page on `OpenKIM <https://openkim.org/>`__ along with
the corresponding BibTex file, and is automatically added to the
DL_POLY_4    *log.cite* file.

Citing the IM software (KIM infrastructure and specific PM or SM codes)
used in the simulation gives credit to the researchers who developed
them and enables open source efforts like OpenKIM to function.