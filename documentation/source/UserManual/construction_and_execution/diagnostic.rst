Warning and Error Processing
============================

The DL_POLY_4 Internal Warning Facility
---------------------------------------

DL_POLY_4 contains a number of various in-built checks scattered
throughout the package which detect a range of possible inconsistencies
or errors. In all cases, such a check fails the subroutine ``warning``
is called, resulting in an appropriate message that identifies the
inconsistency. In some cases an inconsistency is resolved by DL_POLY_4
supplying a default value or assuming a priority of one directive over
the another (in clash of mutually exclusive directives). However, in
other cases this cannot be done and controlled termination of the
program execution is called by the subroutine ``error``. In any case
appropriate diagnostic message is displayed notifying the user of the
nature of the problem.

The DL_POLY_4 Internal Error Facility
-------------------------------------

DL_POLY_4 contains a number of in-built error checks scattered
throughout the package which detect a wide range of possible errors. In
all cases, when an error is detected the subroutine ``error`` is called,
resulting in an appropriate message and termination of the program
execution (either immediately, or after some additional processing). In
some case, if the cause for error is considered to be mendable it is
corrected and the subroutine ``warning`` results in an appropriate
message.

Users intending to insert new error checks should ensure that all error
checks are performed *concurrently* on *all* nodes, and that in
circumstances where a different result may obtain on different nodes, a
call to the global status routine ``gcheck`` is made to set the
appropriate global error flag on all nodes. Only after this is done, a
call to subroutine ``error`` may be made. An example of such a procedure
might be:

    .. code-block::

        > Logical :: safe 
        > safe = (test_condition) 
        > Call gcheck(safe) 
        > If (.not.safe) Call error(message_number)

In this example it is assumed that the logical operation test_condition
will result in the answer *.true.* if it is safe for the program to
proceed, and *.false.* otherwise. The call to ``error`` requires the
user to state the ``message_number`` is an integer which used to
identify the appropriate message to be printed.

.. index:: single: error messages

A full list of the DL_POLY_4 error messages and the appropriate user
action can be found in :ref:`Appendix D<error-messages>`
of this document.
