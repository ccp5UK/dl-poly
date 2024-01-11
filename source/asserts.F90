Module asserts

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! dl_poly_4 module for assertions
  !
  ! copyright - daresbury laboratory
  ! author    - A. Buccheri June 2019
  !
  ! contrib   - J. Wilkins October 2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef SERIAL
  Use mpi,             Only: MPI_COMM_WORLD, MPI_ABORT
#endif
  Use constants, Only : wp
  Use iso_fortran_env, Only: error_unit
  Implicit None

  Integer, Parameter :: error_code_logical = -101
  Real(kind=wp), Parameter :: default_tol = 1.0e-6

  Interface assert
    Module Procedure assert_true
    Module Procedure assert_equal_int
    Module Procedure assert_equal_string
    Module Procedure assert_almost_equal
    Module Procedure assert_almost_equal_rvec
  End Interface assert

  Private
  Public :: assert

Contains

  !> @brief Assert if logical condition is true
  !!
  !! @param[in]     logical_condition    Condition to test
  !! @param[in]     message              Optional message
  !
  Subroutine assert_true(logical_condition, message, passed, passed_accum)
    Logical,          Intent(In   )                 :: logical_condition
    Character(Len=*), Intent(In   ),       Optional :: message
    Logical,          Intent(  Out),       Optional :: passed
    Logical,          Intent(InOut),       Optional :: passed_accum

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

    If (Present(passed)) then
      passed = logical_condition
    End If

    If (Present(passed_accum)) Then
      passed_accum = passed_accum .and. logical_condition
    End If

#endif
  End Subroutine assert_true

  Subroutine assert_equal_int(actual, expected, message, passed, passed_accum)
    Integer,          Intent(In   )           :: actual, expected
    Character(Len=*), Intent(In   ), Optional :: message
    Logical,          Intent(  Out), Optional :: passed
    Logical,          Intent(InOut), Optional :: passed_accum

#ifdef WITH_ASSERT
    If (actual /= expected) Then
      If (Present(message)) Then
        Write (error_unit, '(/,1x,a)') Trim(Adjustl(message))
        Write (error_unit, '(a, i0.1)') "Expected: ", expected
        Write (error_unit, '(a, i0.1)') "Received: ", actual
      End If
!!$#ifdef SERIAL
!!$      Stop error_code_logical
!!$#else
!!$      Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr)
!!$#endif
    End If

    If (Present(passed)) then
      passed = actual == expected
    End If

    If (Present(passed_accum)) Then
      passed_accum = passed_accum .and. (actual == expected)
    End If

#endif
  end Subroutine assert_equal_int

  Subroutine assert_equal_string(actual, expected, message, passed, passed_accum)
    Character(Len=*), Intent(In   )           :: actual, expected
    Character(Len=*), Intent(In   ), Optional :: message
    Logical,          Intent(  Out), Optional :: passed
    Logical,          Intent(InOut), Optional :: passed_accum

#ifdef WITH_ASSERT
    If (actual /= expected) Then
      If (Present(message)) Then
        Write (error_unit, '(/,1x,a)') Trim(Adjustl(message))
        Write (error_unit, '(a, a)') "Expected: ", expected
        Write (error_unit, '(a, a)') "Received: ", actual
      End If
!!$#ifdef SERIAL
!!$      Stop error_code_logical
!!$#else
!!$      Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr)
!!$#endif
    End If

    If (Present(passed)) then
      passed = actual == expected
    End If

    If (Present(passed_accum)) Then
      passed_accum = passed_accum .and. (actual == expected)
    End If

#endif
  end Subroutine assert_equal_string

  Subroutine assert_almost_equal(actual, expected, message, passed, passed_accum, tolerance)
    Real(Kind=wp),    Intent(In   )           :: actual, expected
    Character(Len=*), Intent(In   ), Optional :: message
    Logical,          Intent(  Out), Optional :: passed
    Logical,          Intent(InOut), Optional :: passed_accum
    Real(Kind=wp),    Intent(In   ), Optional :: tolerance

    Real(Kind=wp)                             :: tol

    If (Present(tolerance)) Then
      tol = tolerance
    Else 
      tol = default_tol
    End If


#ifdef WITH_ASSERT
    If (abs(actual - expected) > tol) Then
      If (Present(message)) Then
        Write (error_unit, '(/,1x,a)') Trim(Adjustl(message))
        Write (error_unit, '(a, g15.8)') "Expected: ", expected
        Write (error_unit, '(a, g15.8)') "Received: ", actual
      End If
!!$#ifdef SERIAL
!!$      Stop error_code_logical
!!$#else
!!$      Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr)
!!$#endif
    End If

    If (Present(passed)) then
      passed = abs(actual - expected) < tol
    End If

    If (Present(passed_accum)) Then
      passed_accum = passed_accum .and. (abs(actual - expected) < tol)
    End If

#endif
  end Subroutine assert_almost_equal

  Subroutine assert_almost_equal_rvec(actual, expected, message, passed, passed_accum, tolerance)
    Real(Kind=wp), Dimension(:), Intent(In   )           :: actual, expected
    Character(Len=*),            Intent(In   ), Optional :: message
    Logical,                     Intent(  Out), Optional :: passed
    Logical,                     Intent(InOut), Optional :: passed_accum
    Real(Kind=wp),               Intent(In   ), Optional :: tolerance

    Real(Kind=wp)                                        :: tol

    If (Present(tolerance)) Then
      tol = tolerance
    Else 
      tol = default_tol
    End If

#ifdef WITH_ASSERT
    If (any(abs(actual - expected) > tol)) Then
      If (Present(message)) Then
        Write (error_unit, '(/,1x,a)') Trim(Adjustl(message))
        Write (error_unit, '(a,*(g16.8,","))') "Expected: ", expected
        Write (error_unit, '(a,*(g16.8,","))') "Received: ", actual
        Write (error_unit, '(a,*(g16.8,","))') "Diff:     ", actual-expected
      End If
!!$#ifdef SERIAL
!!$      Stop error_code_logical
!!$#else
!!$      Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr)
!!$#endif
    End If

    If (Present(passed)) then
      passed = All(Abs(actual - expected) < tol)
    End If

    If (Present(passed_accum)) Then
      passed_accum = passed_accum .and. All(Abs(actual - expected) < tol)
    End If

#endif
  end Subroutine assert_almost_equal_rvec

End Module asserts
