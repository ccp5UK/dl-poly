The `execute' Sub-directory

In the supplied version of DL_POLY_4, this sub-directory is intended to
be the `working' directory from which jobs are submitted for execution
and the data files manipulated.  For these purposes the sub-directory
contains a few macros for copying and storing data files from and to
the data sub-directory and for submitting programs for execution.
These are described below.

Executing the Program

To run DL_POLY_4, it is necessary first to ensure that the program is
compiled (from the source sub-directory) and that the files CONTROL,
CONFIG and FIELD are present in the execute subdirectory.  (Some of
the macros described below will help with this.)

In the case of IBM LoadLeveler, parallel execution is performed by
submitting the command:

llsubmit gopoly

where in this example the macro `gopoly', specifies all necessary and
some extra variables controlling the execution.  Execution on one CPU
can simply be performed by typing `./DLPOLY.Z'.  In any case, the
program will allocate the data files automatically.  All output data
files will be returned to the execute sub-directory.

The copy macro

When a job has finished and it needs to be restarted, typing:

copy

in the execute sub-directory will result in some of the output files
being renamed ready for a subsequent run of the program.  Useful if
you can't remember which files to manipulate.

The select macro

If you type, in the execute sub-directory, the command:

select n

where n is a number, the macro will copy the data files needed from
/TESTn/ in the data sub-directory into execute sub-directory.  This is
useful for trying out the test cases stored in the data sub-directory.

The store macro

If you type, in the execute sub-directory, the command:

store n

where n is an alphanumeric string, the macro will store all input and
output files within the execute directory into the data sub-directory
as the corresponding files /TESTn/.

The cleanup macro

After a job has finished, you can clean up the execute sub-directory
by typing:

cleanup

in the sub-directory.  This will erase all unessential files prior to
the next run of the program.  Do not use this macro before you have
stored all your useful data!
