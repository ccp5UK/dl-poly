Module unit_test

  Use comms,              Only: comms_type
  Use test_configuration, Only: run_configuration_tests
  Use test_units,         Only: run_units_tests
  Use test_control,       Only: run_control_tests
  Use test_vdw,           Only: run_vdw_tests
  Use test_integrators,   Only: run_integrators_tests
  Use test_hash,          Only: run_hash_tests
  Implicit None

  !> Logicals indicating whether tests should be run for
  !> corresponding module. Add as required
  Type, Public :: testing_type
    Logical, Public :: configuration = .false.
    Logical, Public :: units = .false.
    Logical, Public :: control = .false.
    Logical, Public :: dftb_library = .false.
    Logical, Public :: vdw = .false.
    Logical, Public :: integrators = .false.
    Logical, Public :: hash = .false.
  Contains
    Procedure :: all => set_all_tests_true
    Procedure :: run => run_unit_tests
  End Type testing_type

Contains

  Subroutine set_all_tests_true(this)
    Class(testing_type), Intent(inout) :: this
    this%configuration = .true.
    this%units = .true.
    this%control = .true.
    this%dftb_library = .true.
    this%vdw = .true.
    this%integrators = .true.
    this%hash = .true.
  End Subroutine set_all_tests_true

  Subroutine run_unit_tests(this, comm, passed_all)
    Class(testing_type), Intent(InOut) :: this
    Type(comms_type)                   :: comm
    Logical,             Intent(  Out) :: passed_all
    Logical                            :: passed = .true.

    passed_all = .true.

    If (this%units) Then 
      Call run_units_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%control) Then 
      Call run_control_tests(comm, passed)
      passed_all = passed_all .and. passed
    End If
    
    If (this%configuration) Then
      Call run_configuration_tests(comm, passed)
      passed_all = passed_all .and. passed
    End If
    
    If (this%vdw) Then 
      Call run_vdw_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%integrators) Then
      Call run_integrators_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%hash) Then
      Call run_hash_tests(passed)
      passed_all = passed_all .and. passed
    End If

  End Subroutine run_unit_tests
end Module unit_test
