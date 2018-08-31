Release notes for DL_POLY_4
============================

Version 4.09 August 2018
-------------------------

+ new potentials
  + modified morse 
  + Ziegler-Biersack-Littmark, ZBL, ZBL mixed with Morse and Buckingham
  + Mei-Devenport-Fernando taper for Lennard-Jones, 12-6 and  Buckingham 
  + Rydberg
  + modified LJ
+ rdf now can calculate error bars
+ nfold deals correctly with cubic cells
+ CONTROL file can be passed as command line argument under any names
+ main input/output filenames can be changed via directives from control file
+ fix restart for velocity autocorrelation functions
+ fix reading from tables for bonds,angles and dihedrals, thanks to Tom Potter from University of Durham
+ two temperature model
+ fix bugs in trajectory writing
+ doubled number of tests
+ windows crosscompilation from linux both serial and parallel
+ fix bug in units for TABLE file
+ enhance STATIS writing thanks to Mike Pacey, University of Lancaster
