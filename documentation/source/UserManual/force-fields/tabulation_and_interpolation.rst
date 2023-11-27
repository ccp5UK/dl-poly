Tabulation and interpolation in the treatment of intermolecular interactions
============================================================================

By default DL_POLY_4 tabulates in memory most of the intermolecular
interactions keeping values of the potential and the negative of its
first derivative times the distance (or virial) over an equidistant
grid. This is done for reasons of speed as due to the large variety of
potential forms, some could be quite expensive to evaluate if run
unoptimised. The memory tabulation could be overridden for non-tabulated
interactions upon user specified options such as **metal direct** for
metal interactions and **vdw direct** for van der Waals interactions.

When energy and force are calculated for tabulated interactions a
3-point interpolation scheme of our own is used to interpolate the value
for the requested distance.

A 5-point interpolation is used for finding the numerical derivatives of
(2B)(E)EAM type potentials which are supplied in TABEAM by the user. For
this a Lagrange formula is used, which can be found in any textbook on
numerical methods.