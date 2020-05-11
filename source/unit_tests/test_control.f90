Program test_control

  use comms, only : comms_type
  use filename, only : file_type
  use new_control, only : read_new_control
  Use errors_warnings, only : error
  use control_parameter_module, only : parameters_hash_table, print_set
  use units, only : initialise_units, set_timestep, convert_units

  type( parameters_hash_table ) :: params
  type(comms_type) :: comm
  type(file_type) :: control_file
  Logical :: can_parse
  Real(kind=8) :: ts


  control_file = file_type(filename='test_new_control')

  call initialise_units()
  call read_new_control(control_file, params, comm, can_parse)
  if (.not. can_parse) call error(0, 'Cannot parse control')
  call params%retrieve('timestep', ts)
  call set_timestep(ts)

  call print_set(params)
!  call params%help()

end Program test_control
