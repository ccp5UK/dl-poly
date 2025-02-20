Module unit_test

  Use comms,              Only: comms_type
  Use test_configuration, Only: run_configuration_tests
  Use test_smearing,      Only: run_smearing_test
  Use test_units,         Only: run_units_tests
  Use test_control,       Only: run_control_tests
  Use test_vdw,           Only: run_vdw_tests
  Use test_integrators,   Only: run_integrators_tests
  Use test_hash,          Only: run_hash_tests
  Use test_parse,         Only: run_parse_tests
  Use test_numerics,      Only: run_numerics_tests
  Use test_mdpd_ld,       Only : run_mdpd_tests
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
    Logical, Public :: parse = .false.
    Logical, Public :: numerics = .false.
    Logical, Public :: mdpd = .false.
    Logical, Public :: smearing = .false.
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
    this%parse = .true.
    this%numerics = .true.
    this%mdpd = .true.
    this%smearing = .true.
  End Subroutine set_all_tests_true

  Subroutine run_unit_tests(this, comm, eu)
    Class(testing_type), Intent(InOut) :: this
    Type(comms_type),    Intent(InOut) :: comm
    Integer,             Intent(In   ) :: eu

    Logical :: passed_all
    Logical :: passed = .true.

    passed_all = .true.

    If (this%units) Then
      Write(eu, '(a)') "Running test: units"
      Call run_units_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%control) Then 
      Write(eu, '(a)') "Running test: control"
      Call run_control_tests(comm, passed)
      passed_all = passed_all .and. passed
    End If
    
    If (this%configuration) Then
      Write(eu, '(a)') "Running test: configuration"
      Call run_configuration_tests(comm, passed)
      passed_all = passed_all .and. passed
    End If
    
    If (this%vdw) Then 
      Write(eu, '(a)') "Running test: vdw"
      Call run_vdw_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%integrators) Then
      Write(eu, '(a)') "Running test: integrators"
      Call run_integrators_tests(passed)
      passed_all = passed_all .and. passed
    End If
    
    If (this%mdpd) Then 
      Write(eu, '(a)') "Running test: mdpd"
      Call run_mdpd_tests(passed)
      passed_all = passed_all .and. passed 
    End If

    If (this%smearing) Then 
      Write(eu, '(a)') "Running test: charge smearing"
      Call run_smearing_test(passed)
      passed_all = passed_all .and. passed
    End If  

    If (this%hash) Then
      Write(eu, '(a)') "Running test: hash"
      Call run_hash_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%parse) Then
      Write(eu, '(a)') "Running test: parse"
      Call run_parse_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (this%numerics) Then
      Write(eu, '(a)') "Running test: numerics"
      Call run_numerics_tests(passed)
      passed_all = passed_all .and. passed
    End If

    If (passed_all) Then
      Write(eu, '(a)') "Status: PASSED"
    Else
      Write(eu, '(a)') "Status: FAILED"
    End If

  End Subroutine run_unit_tests
end Module unit_test
