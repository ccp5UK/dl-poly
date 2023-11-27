.. _macros:

Appendix C: DL_POLY_4 Macros 
++++++++++++++++++++++++++++

.. index:: single: sub-directory

Macros are simple executable files containing standard UNIX commands. A
number of the are supplied with DL_POLY_4 and are found in the *execute*
sub-directory. These are not guaranteed to be immaculate but with little
adaptation they can become a useful tool to a researcher. The available
macros are as follows:

-  cleanup

-  copy

-  gopoly

-  gui

-  select

-  store

The function of each of these is described below. It is worth noting
that most of these functions could be performed by the DL_POLY Java GUI
:cite:`smith-gui`.

*cleanup*
~~~~~~~~~

cleanup removes several standard data files from the execute
:index:`sub-directory`. It contains the UNIX commands:

::

   rm -vf OUTPUT STATIS REVCON RVEIVE CFGMIN DEFECTS *DAT* *PMF *TAB MSDTMP

It is useful for cleaning the :index:`sub-directory` up after a run. (Useful data
should be stored elsewhere however!)

*copy*
~~~~~~

copy invokes the UNIX commands:

::

   mv -v CONFIG CONFIG.OLD
   mv -v REVCON CONFIG
   mv -v REVIVE REVOLD

which collectively prepare the DL_POLY_4 files in the *execute*
:index:`sub-directory` for the continuation of a simulation. It is always a good
idea to store these files elsewhere in addition to using this macro.

*gopoly*
~~~~~~~~

gopoly is used to submit a DL_POLY_4 job to the HPC\ :math:`x`, which
operates a LOAD-LEVELER job queuing system. It invokes the following
script:

::

   #@ shell = /usr/bin/tcsh
   #
   #@ job_type = parallel
   #@ job_name = gopoly
   #
   #@ cpus = 32
   #
   #@ node_usage = not_shared
   #@ network.MPI = csss,shared,US
   #
   #@ wall_clock_limit = 00:30:00
   #@ account_no = mine
   #
   #@ output = $(job_name).$(schedd_host).$(jobid).out
   #@ error  = $(job_name).$(schedd_host).$(jobid).err
   #@ notification = never
   #
   #@ bulkxfer = yes
   #@ data_limit = 850000000
   #@ stack_limit = 10000000
   #
   #@ queue
   #
   # ENVIRONMENT SETTINGS
   #
   setenv MP_EAGER_LIMIT 65536
   setenv MP_SHARED_MEMORY yes
   setenv MEMORY_AFFINITY MCM
   setenv MP_TASK_AFFINITY MCM
   setenv MP_SINGLE_THREAD yes
   #
   poe  ./DLPOLY.Z

Using LOADLEVELLER, the job is submitted by the UNIX command:
llsubmit gopoly
where llsubmit is a local command for submission to the IBM SP4
cluster. The number of required nodes and the job time are indicated
in the above script.

*gui*
~~~~~

gui is a macro that starts up the DL_POLY_4 Java GUI. It invokes the
following UNIX commands:

::

   java -jar ../java/GUI.jar $1 &

In other words the macro invokes the Java Virtual Machine which executes
the instructions in the Java archive file GUI.jar, which is stored in
the *java* subdirectory of . (Note: Java 1.3.0 or a higher version is
required to run the GUI.)

*select*
~~~~~~~~

select is a macro enabling easy selection of one of the test cases. It
invokes the UNIX commands:

::

   cp -vpLH ../data/TEST$1/CONTROL   .
   cp -vpLH ../data/TEST$1/CONFIG    .
   cp -vpLH ../data/TEST$1/HISTORY   .
   cp -vpLH ../data/TEST$1/FIELD     .
   cp -vpLH ../data/TEST$1/MPOLES    .
   cp -vpLH ../data/TEST$1/TAB*      .
   cp -vpLH ../data/TEST$1/REFERENCE .
   cp -vpLH ../data/TEST$1/Ce.dat    .
   cp -vpLH ../data/TEST$1/g.dat     .

requires one argument (an integer) to be specified:
select n
where n is test case number, which ranges from 1 to 18.

This macro sets up the required input files in the *execute*
:index:`sub-directory` to run the n-th test case. The last three copy commands
may not be necessary in most cases.

*store*
~~~~~~~

The store macro provides a convenient way of moving data back from the
*execute* :index:`sub-directory` to the *data* :index:`sub-directory`. It invokes the UNIX
commands:

::

   mkdir -pv          ../data/TEST$1
   cp -vpLH CONTROL   ../data/TEST$1
   cp -vpLH CONFIG    ../data/TEST$1
   cp -vpLH FIELD     ../data/TEST$1
   cp -vpLH MPOLES    ../data/TEST$1
   cp -vpLH TAB*      ../data/TEST$1
   cp -vpLH REFERENCE ../data/TEST$1
   cp -vpLH HISTORY   ../data/TEST$1
   cp -vpLH Ce.dat    ../data/TEST$1
   cp -vpLH g.dat     ../data/TEST$1
   mv -v     OUTPUT   ../data/TEST$1
   mv -v     STATIS   ../data/TEST$1
   mv -v     REV*     ../data/TEST$1
   mv -v     CFGMIN   ../data/TEST$1
   mv -v     HISTORF  ../data/TEST$1
   mv -v     DEFECTS  ../data/TEST$1
   mv -v     *DAT*    ../data/TEST$1
   mv -v     *PMF     ../data/TEST$1
   mv -v     *TAB     ../data/TEST$1
   mv -v     MSDTMP   ../data/TEST$1
   mv -v     DUMP_E   ../data/TEST$1
   mv -v     LATS_*   ../data/TEST$1
   mv -v     PEAK_*   ../data/TEST$1
   chmod -R a-w       ../data/TEST$1

which first creates a new DL_POLY *data/TEST..* :index:`sub-directory` and then
moves the standard output data files into it.

store requires one argument:
``store n``
where ``n`` is a unique string or number to label the output data in the
*data/TESTn* :index:`sub-directory`.

Note that store sets the file access to read-only. This is to prevent
the store macro overwriting existing data without your knowledge.
