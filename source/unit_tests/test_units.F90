Module test_units

  Use asserts, Only: assert
  Use kinds,   Only: STR_LEN,&
                     wp
  Use units,   Only: convert_units,&
                     destroy_units,&
                     init_scheme,&
                     initialise_units,&
                     set_out_units,&
                     to_out_units,&
                     units_scheme

  Implicit None

Contains

  Subroutine run_units_tests(passed)
    Logical, Intent(  Out)   :: passed
    Real(kind=wp), Parameter :: atm_GPa = 0.00010132501_wp, Cal_J = 4.1842_wp, &
                                eV_J = 1.602176634e-19_wp, Ha_eV = 27.211386245988_wp, &
                                Ha_J = 4.3597482e-18_wp, kB_JK = 1.3804703422756560E-023, &
                                pound_kg = 0.45359237_wp

    Character(Len=STR_LEN) :: unit_name
    Real(kind=wp)          :: unit_test
    Type(units_scheme)     :: scheme

!J
!eV
    passed = .true.

    Call initialise_units()

    unit_test = convert_units(1.0_wp, 'e_h', 'e.V')
    Call assert(unit_test, Ha_eV, "Hartree->eV fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'cal', 'J')
    Call assert(unit_test, Cal_J, "Cal->J fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'e_h', 'J')
    Call assert(unit_test, Ha_J, "Hartree->J fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'e.V', 'J')
    Call assert(unit_test, eV_J, "eV->J fail", passed_accum = passed)
    
    unit_test = convert_units(10.0_wp, 'J/mol', 'internal_e')
    Call assert(unit_test, 1.0_wp, "J/mol->Internal fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'V', 'J/C')
    Call assert(unit_test, 1.0_wp, "V->J/C fail", passed_accum = passed)   

    unit_test = convert_units(1.0_wp, 'ang/ps', 'internal_v')
    Call assert(unit_test, 1.0_wp, "ang/ps->internal_v fail", passed_accum = passed)   

    unit_test = convert_units(1.0_wp, 'atm', 'GPa')
    Call assert(unit_test, atm_GPa, "atm->GPa fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'N/m^2', 'Pa')
    Call assert(unit_test, 1.0_wp, "N/m^2->Pa fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'grav', 'm/s^2')
    Call assert(unit_test, 9.81_wp, "g->m/s^2 fail", passed_accum = passed)

    unit_test = convert_units(1.0_wp, 'J', 'kg.m^2/s^2')
    Call assert(abs(unit_test-1.0_wp) < 1e-3_wp, "J->kg.m^2/s^2 fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'J', 'W.s')
    Call assert(unit_test, 1.0_wp, "J->Ws fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'J', 'N.m')
    Call assert(abs(unit_test-1.0_wp) < 1e-3_wp, "J->N.m fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'lb', 'kg')
    Call assert(unit_test, pound_kg, "lb->kg fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'ft', 'in')
    Call assert(unit_test, 12.0_wp, "ft->in fail", passed_accum = passed)
    
    unit_test = convert_units(1.0_wp, 'k_B', 'm^2.kg.s^-2/K')
    Call assert(unit_test, kB_JK, "k_b->m^2.kg.s^-2/K fail", passed_accum = passed)

    Call to_out_units(1.0_wp, 'Ang', unit_test, out_unit=unit_name)
    Call assert(unit_test, 1.0_wp, "Out units internal fail", passed_accum = passed)
    
    Call assert(unit_name, 'internal_l', "Out units internal fail", passed_accum = passed)
    

    scheme = init_scheme( &
             length='m', &
             time='s', &
             mass='kg', &
             charge='C', &
             energy='J', &
             pressure='Pa', &
             force='N', &
             velocity='m/s', &
             power='W', &
             surf_ten='N/m', &
             emf='V')
    Call set_out_units(scheme)

    Call to_out_units(1.0_wp, 'Ang', unit_test, out_unit=unit_name)
    Call assert(unit_test, 1.0e-10_wp, "Out units internal fail", passed_accum = passed)
    
    Call assert(unit_name, 'm', "Out units internal fail", passed_accum = passed)
    
    Call destroy_units()

  End Subroutine run_units_tests

End Module test_units
