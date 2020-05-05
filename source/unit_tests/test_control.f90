program test_control

  use comms, only : comms_type
  use new_control, only : initialise_control, read_new_control => parse_file
  use control_parameter_module, only : parameters_hash_table, print_set
  use units, only : initialise_units, set_timestep, convert_units

  type( parameters_hash_table ) :: params
  type(comms_type) :: comm
  Real(kind=8) :: ts
  Character(Len = 50) :: line
  Integer :: test

  call initialise_units()
  open(newunit=test, file='test_new_control', status='old', action='read')
  call initialise_control(params)
  call read_new_control(test, params, comm)

  call params%retrieve('timestep', ts)
  call set_timestep(ts)
  call print_set(params)
!  call params%help()

end program test_control
