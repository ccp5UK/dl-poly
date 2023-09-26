Module unit_test

  Use comms, only : comms_type
  Use test_configuration, only : run_configuration_tests
  Use test_units, only : run_units_tests
  Use test_control, only : run_control_tests
  Use test_vdw, Only : run_vdw_tests
  Use test_integrators, Only : run_integrators_tests
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
  End Subroutine set_all_tests_true

  Subroutine run_unit_tests(this, comm)
    Class(testing_type), Intent(inout) :: this
    Type(comms_type) :: comm

    if (this%units) Call run_units_tests()
    if (this%control) Call run_control_tests(comm)
    if (this%configuration) Call run_configuration_tests(comm)
    if (this%vdw) Call run_vdw_tests()
    if (this%integrators) Call run_integrators_tests()
  End Subroutine run_unit_tests
end Module unit_test
