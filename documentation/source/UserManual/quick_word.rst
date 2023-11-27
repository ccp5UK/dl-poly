Quick Word / INSTALL & RUN
==========================

**For the experienced and quick minded this is a very brief resume
of how to INSTALL & RUN DL_POLY_4 (which is no excuse for skipping the
Introduction, Chapter** :ref:`introduction_ch`\ **!).  For the rest of us it
sketches out how to start running DL_POLY_4 jobs and where one should look
to obtain more detailed information if need be.**

If you have followed the procedure for obtaining and downloading the
DL_POLY_4 package (see Obtaining~the~Source~Code, Section~:ref:`source-code-sec`),
have successfully unpacked it and are ready to compile the source code,
then jump to the INSTALL Notes in the ``INSTALL`` file, both in the main
distribution directory as well as in :ref:`Appendix E<readme>`.

If you have compiled successfully then a freshly date-stamped file,
named DLPOLY.Z 
should appear in the listing of the *execute* subdirectory (the
'``ls -haltr``' command issued on a Linux/Unix-like shell within
*execute* will place the executable in the last row of the list).
If **unsuccessful** then you should read the Section: 
:ref:`compilation`

To run the code you first need to place the necessary input files
within *execute*.  TEST cases containing suitable input files,
as well as examples of output files, can be obtained by running CMake 
or building DL_POLY_4 with the CMake option ``BUILD_TESTING=ON``. For example

    .. code-block::
        
        folder="build-mpi-testing"
        rm -rf $folder && mkdir $folder && pushd $folder
        cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON 

Examine the contents of *data*/``README.txt`` and
*bench*/``README.txt`` for more information.  To run the serial
version you simply issue the command ``DLPOLY.Z``
within the *execute* subdirectory.  If you have compiled a parallel
version and are running it on a parallel machine, then naturally you
will need to familiarise yourself with the local procedures of how to
run jobs on that machine.  In general though, running a parallel job
will usually require that you issue a necessary run command (e.g.
``mpirun -n 8 DLPOLY.Z``) or submit a job script from within *execute*.

If you need to know more then search the manual and examine sections of
interests.  You may also wish to visit DL\_POLY project web-page `<http://www.ccp5.ac.uk/DL\_POLY/>`_
and examine the useful links within (FAQ, User Forum, etc.).

If you are looking to gain more in depth experience, then regular training
workshops are available.  To find about upcoming workshops, subscribe to
our mail list by following instructions in Section :ref:`otherInfoSection`.

If you need one-to-one training, wish to collaborate scientifically and/or
would like to become a contributor/developer then get in touch with me,
Dr. I.T. Todorov, by emailing to
`ilian.todorov@stfc.ac.uk<mailto:ilian.todorov@stfc.ac.uk>`_.

Best of luck!