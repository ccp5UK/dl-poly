program test_control

  use comms, only : comms_type
  use new_control, only : parameters_hash_table, initialise_control, print_set, read_new_control => parse_file
  use units, only : initialise_units, set_timestep, convert_units

  type( parameters_hash_table ) :: params
  type(comms_type) :: comm
  Real(kind=8) :: ts
  Integer :: test


  call initialise_units()
  ts = convert_units(1.0_8, 'fs', 's')
  print*, ts
  open(newunit=test, file='test_new_control')
  call initialise_control(params)
  call read_new_control(test, params, comm)
  call params%retrieve('timestep', ts)
  call set_timestep(ts)
  call print_set(params)

end program test_control
