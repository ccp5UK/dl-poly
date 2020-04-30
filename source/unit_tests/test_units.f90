program test_units

  Use asserts, only : assert
  Use constants, only : tenunt, boltz
  Use units, only : initialise_units, convert_units
  Use kinds, only : wp
  Implicit None

  Real( kind = wp ), parameter :: eV_J = 1.602176634e-19_wp !J
  Real( kind = wp ), parameter :: Ha_eV = 27.211386245988_wp !eV

  Real( kind = wp ) :: unit_test

  call initialise_units()

  unit_test = convert_units(1.0_wp, 'e_h', 'e.V')
  print*, "Ha -> eV", unit_test

  unit_test = convert_units(1.0_wp, 'cal', 'J')
  print*, "Cal -> J", unit_test

  unit_test = convert_units(1.0_wp, 'e_h', 'J')
  print*, "Ha -> J", unit_test

  unit_test = convert_units(1.0_wp, 'e.V', 'J')
  print*, "eV -> J", unit_test

  unit_test = convert_units(10.0_wp, 'J/mol', 'internal_e')
  print*, "10 J/mol -> Internal", unit_test

  unit_test = convert_units(1.0_wp, 'V', 'J/C')
  print*, "V -> J/C", unit_test

  unit_test = convert_units(1.0_wp, 'm/s', 'internal_v')
  print*, "m/s -> ang/ps", unit_test

  unit_test = convert_units(1.0_wp, 'atm', 'Gpa')
  print*, "atm -> GPa", unit_test

  print*, 1.0e-12_wp * 1.0e-7_wp / (boltz * tenunt) / 3.0_wp
  unit_test = convert_units(1.0_wp / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')
  print*, unit_test

end program test_units
