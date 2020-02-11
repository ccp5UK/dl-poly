Module asserts

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! dl_poly_4 module for assertions
  !
  ! copyright - daresbury laboratory
  ! author    - A. Buccheri June 2019
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef SERIAL
  Use mpi,             Only: MPI_COMM_WORLD, MPI_ABORT
#endif
  Use iso_fortran_env, Only: error_unit
  Implicit None

  Integer, Parameter :: error_code_logical = -101

  Interface assert
    Module Procedure assert_true
  End Interface assert

  Private
  Public :: assert

Contains

  !> @brief Assert if logical condition is true
  !!
  !! @param[in]     logical_condition    Condition to test
  !! @param[in]     message              Optional message
  !
  Subroutine assert_true(logical_condition, message)
    Logical, Intent(In)           :: logical_condition
    Character(Len=*), Intent(In), Optional :: message
    Integer                                    :: ierr
#ifdef WITH_ASSERT
    If (.not. logical_condition) Then
      If (Present(message)) Then
        Write (error_unit, '(/,1x,a)') Trim(Adjustl(message))
      End If
!!$#ifdef SERIAL
!!$      Stop error_code_logical
!!$#else
!!$      Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr)
!!$#endif
    End If

#endif
  End Subroutine assert_true

End Module asserts
