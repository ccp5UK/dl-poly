program test_control

  use new_control, only : parameters_hash_table, initialise_control

  type( parameters_hash_table ) :: params

  call initialise_control(params)
  call params%help('ewald_nsplines')


end program test_control
