DL_POLY_4 Java GUI
==================

The DL_POLY_4 adapted GUI is provided on "AS IS" terms under the
Licence (within this source directory), i.e. with no support what-
so-ever!

Shell commands within this source directory:
============================================
build - use this to rebuild the GUI on your computer, but only of you
        have the Java SDK installed on your computer
try   - script to run the program from within this directory, even if
        the script is executed as a soft link anywhere on the OS
rmc   - script to remove the class files, once the GUI.jar file has
        been built (by using the "build" script for example)

Note that the GUI.jar file is supposed to be universally exportable -
it does not have to be inside the DL_POLY_4 directory to have most of
the functionality working.  However, to run any locally stored test
cases, it is necessary to be in the default DL_POLY_4 directory
structure.  The GUI invokes DLPOLY.Z directly (see the file
Execute.java), but one could in principle submit a script instead,
to send the job to remote machines or run in parallel locally.  In
the latter casess one will lose the ability to find out what the job
status is, but it would still submit the job OK.

Good luck!

Ilian Todorov & Bill Smith
