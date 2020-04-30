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

  unit_test = convert_units(1.0_wp, 'N/m^2', 'Pa')
  print*, "N/m^2 -> Pa", unit_test

  unit_test = convert_units(1.0_wp, 'grav', 'm/s^2')
  print*, "g -> m/s^2", unit_test

  unit_test = convert_units(1.0_wp, 'J', 'kg.m^2/s^2')
  print*, 'J -> kgm^2/s^2', unit_test

  unit_test = convert_units(1.0_wp, 'J', 'W.s')
  print*, 'J -> Ws', unit_test

  unit_test = convert_units(1.0_wp, 'J', 'N.m')
  print*, 'J -> Nm', unit_test

  unit_test = convert_units(1.0_wp, 'lb', 'kg')
  print*, 'lb -> kg', unit_test

  unit_test = convert_units(1.0_wp, 'ft', 'in')
  print*, 'ft -> "', unit_test

  unit_test = convert_units(1.0_wp, 'k_B', 'm^2.kg.s^-2/K')
  print*, 'K_B in J/K', unit_test

  print*, 1.0e-7_wp / (boltz * tenunt)
  unit_test = convert_units(1.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')
  print*, unit_test

end program test_units
