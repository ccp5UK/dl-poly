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
  Real(kind=wp), Parameter :: tol = 1.0e-6

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
  Subroutine assert_true(logical_condition, message, passed)
    Logical, Intent(In)           :: logical_condition
    Character(Len=*), Intent(In), Optional :: message
    Logical, Intent(Out), Optional :: passed

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

    if (present(passed)) then
      passed = .not. logical_condition
    end if

#endif
  End Subroutine assert_true

  Subroutine assert_equal_int(actual, expected, message, passed)
    Integer, Intent(In) :: actual, expected
    Character(Len=*), Intent(In), Optional :: message
    Logical, Intent(Out), Optional :: passed

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

    if (present(passed)) then
      passed = actual == expected
    end if
#endif
  end Subroutine assert_equal_int

  Subroutine assert_equal_string(actual, expected, message, passed)
    Character(Len=*), Intent(In) :: actual, expected
    Character(Len=*), Intent(In), Optional :: message
    Logical, Intent(Out), Optional :: passed

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

    if (present(passed)) then
      passed = actual == expected
    end if
#endif
  end Subroutine assert_equal_string

  Subroutine assert_almost_equal(actual, expected, message, passed)
    Real(Kind=wp), Intent(In) :: actual, expected
    Character(Len=*), Intent(In), Optional :: message
    Logical, Intent(Out), Optional :: passed

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

    if (present(passed)) then
      passed = abs(actual - expected) < tol
    end if
#endif
  end Subroutine assert_almost_equal

  Subroutine assert_almost_equal_rvec(actual, expected, message, passed)
    Real(Kind=wp), Dimension(:), Intent(In) :: actual, expected
    Character(Len=*), Intent(In), Optional :: message
    Logical, Intent(Out), Optional :: passed

#ifdef WITH_ASSERT
    If (any(abs(actual - expected) > tol)) Then
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

    if (present(passed)) then
      passed = All(abs(actual - expected) < tol)
    end if
#endif
  end Subroutine assert_almost_equal_rvec

End Module asserts
