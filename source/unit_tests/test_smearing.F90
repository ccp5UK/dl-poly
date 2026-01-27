Module test_smearing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for unit testing the charge smearing functions. Tests
  ! standalone smearing correction functions as well as the generated erfc
  ! tables used for Coulombic interactions. No smearing case is equivalent to
  ! a pure real space SPME based coulomb interaction contribution.
  !
  ! copyright - daresbury laboratory
  ! author    - b.t.speake July 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use asserts,            Only: assert
  Use charge_smearing,    Only: linear_smearing, smearing_correction, slater_exp_smearing, &
                                slater_apprx_smearing
  Use kinds,              Only: wp
  Use electrostatic,      Only: electrostatic_type, ELECTROSTATIC_SPME, SMEARING_LINEAR, &
                                SMEARING_SLATER_EXP, SMEARING_SLATER_TRUNCATED, &
                                SMEARING_GAUSSIAN


  Contains

  Subroutine run_smearing_test(passed)
    Logical, Intent(InOut) :: passed
    passed = .true.
    Call no_smearing_test(passed)
    Call linear_smearing_test(passed)
    Call slater_smearing_test(passed)
    Call slater_apprx_smearing_test(passed)
    Call gaussian_smearing_test(passed)

  End Subroutine run_smearing_test

  Subroutine no_smearing_test(passed)
    Logical, Intent(InOut) :: passed
    Type(electrostatic_type) :: electro
    Real(Kind=wp), dimension(10) :: pot_tab, frc_tab

    pot_tab = [&
      1.105325582682910E-002_wp,  7.293980936342976E-007_wp,  3.095634086667356E-013_wp, &
      6.369616962263265E-022_wp,  5.755815839106417E-033_wp,  2.180462812244352E-046_wp, &
      3.376009513850959E-062_wp,  2.103945172770762E-080_wp,  5.225992287472692E-101_wp, &
      5.139379529340524E-124_wp &
    ]
    frc_tab = [&
      2.923630207768757E-002_wp,  1.583406605716887E-006_wp,  6.392171975873630E-013_wp, &
      1.280993235932933E-021_wp,  1.134133185413451E-032_wp,  4.219460859476861E-046_wp, &
      6.425938515205286E-062_wp,  3.944343392416672E-080_wp,  9.661830082713158E-101_wp, &
      9.381209698089542E-124_wp &
    ]

    electro%key = ELECTROSTATIC_SPME
    Call electro%init_erf_tables(10)
    Call electro%erfcgen(10.0_wp, 1.0_wp)

    Call assert(electro%erfc_over_r%table, pot_tab, "No smearing (pure spme) energies differ from expected", passed_accum=passed)
    Call assert(electro%erfc_gamma%table, frc_tab, "No smearing (pure spme) force contributions differ from expected", &
                  passed_accum=passed)

  End Subroutine no_smearing_test

  Subroutine linear_smearing_test(passed)
    Type(electrostatic_type) :: electro
    Type(smearing_correction) :: smear
    Real(Kind=wp) :: r_ij, r_sm, r_cut
    Real(Kind=wp), Dimension(10) :: pot_tab, frc_tab
    Logical, Intent(InOut) :: passed

    r_ij = 2.0_wp
    r_sm = 2.5_wp
    r_cut = 6.0_wp

    pot_tab = [-0.200823686687457_wp,       0.253130568583740_wp,      -5.251976538909083E-003_wp, &
               -8.217339933034941E-005_wp,  3.095634086667356E-013_wp,  3.630860707313369E-018_wp, &
                6.091771272421279E-024_wp,  1.440933202113285E-030_wp,  4.759289549312803E-038_wp, &
                2.180462812244352E-046_wp]
    frc_tab = [-0.311316293071378_wp,       2.651005403186380E-002_wp, -5.243053583648856E-003_wp, &
               -1.232577170563967E-004_wp,  6.392171975873630E-013_wp,  7.371138462006662E-018_wp, &
                1.219771565454547E-023_wp,  2.850181451370495E-030_wp,  9.308230029307534E-038_wp, &
                4.219460859476861E-046_wp]

    electro%key = ELECTROSTATIC_SPME
    electro%smear = SMEARING_LINEAR
    electro%r_smear = r_sm
    Call electro%init_erf_tables(10)
    Call electro%erfcgen(r_cut, 1.0_wp)

    Call assert(electro%erfc_over_r%table, pot_tab, "Linear smearing energies differ from expected", passed_accum=passed)
    Call assert(electro%erfc_gamma%table, frc_tab, "Linear smearing forces contribs differ from expected", &
                  passed_accum=passed)


    smear = linear_smearing(r_ij, r_sm)

    Call assert(smear%energy,-0.501583276748699_wp, "Linear smearing energy (direct) differs from expected", &
                  passed_accum=passed)
    Call assert(smear%force, -0.166068601127945_wp, "Linear smearing force contrib (direct) differs from expected", &
                  passed_accum=passed)

    r_sm = 1.5_wp
    smear = linear_smearing(r_ij, r_sm)
    Call assert(smear%energy, 5.871146576667304E-003_wp, "Linear smearing energy (direct) differs from expected R=1.5", &
                  passed_accum=passed)
    Call assert(smear%force, 7.045375801336790E-002_wp, "Linear smearing force contrib (direct) differs from expected R=15", &
                  passed_accum=passed)

    r_ij = 5.0_wp
    smear = linear_smearing(r_ij, r_sm)
    Call assert(smear%energy, 0.0_wp, "Linear smearing energy (direct) differs from expected r=5", &
                  passed_accum=passed)
    Call assert(smear%force, -0.0_wp, "Linear smearing force contrib (direct) differs from expected r=5", &
                  passed_accum=passed)
  End Subroutine

  Subroutine slater_smearing_test(passed)
    Logical, Intent(InOut) :: passed
    Type(electrostatic_type) :: electro
    Type(smearing_correction) :: smear

    Real(Kind=wp), Dimension(10) :: pot_tab, frc_tab

    pot_tab = [ &
    -0.302078009803406_wp,      -6.578145370473336E-002_wp, -1.580610066837843E-002_wp, &
    -3.834987691167041E-003_wp, -9.174324247674806E-004_wp, -2.154604999977927E-004_wp, &
    -4.970345488738767E-005_wp, -1.128244117311474E-005_wp, -2.524968737142154E-006_wp, &
    -5.581037581401504E-007_wp &
    ]
    frc_tab = [ &
    -0.170321265409243_wp,      -1.720319082936640E-002_wp, -2.681702904996800E-003_wp, &
    -4.907240587027374E-004_wp, -9.508354135709846E-005_wp, -1.884853750639163E-005_wp, &
    -3.770595444284361E-006_wp, -7.566075758087742E-007_wp, -1.518521778402211E-007_wp, &
    -3.044131350035167E-008_wp &
    ]

    electro%key = ELECTROSTATIC_SPME
    electro%smear = SMEARING_SLATER_EXP
    electro%r_smear = 2.0_wp
    electro%b_smear = 0

    Call electro%init_erf_tables(10)
    Call electro%erfcgen(10.0_wp, 1.0_wp)

    Call assert(electro%erfc_over_r%table, pot_tab, "Slater (exp) smearing tabulated energies differ from expected", &
                passed_accum=passed)
    Call assert(electro%erfc_gamma%table, frc_tab,"Slater (exp) smearing tabulated forces differ from expected", &
                passed_accum=passed)

    smear = slater_exp_smearing(4.0_wp, 2.0_wp, 0)
    Call assert(smear%energy, 0.148051414350601_wp, &
              "Slater (exp) smearing direct energy differs from expected", passed_accum=passed)
    Call assert(smear%force, 0.506732675921646_wp, &
              "Slater (exp) smearing direct force differs from expected", passed_accum=passed)

  End Subroutine

  Subroutine slater_apprx_smearing_test(passed)
    Logical, Intent(InOut) :: passed
    Type(electrostatic_type) :: electro
    Type(smearing_correction) :: smear

    Real(Kind=wp), Dimension(10) :: pot_tab, frc_tab

    pot_tab = [&
    -0.196709907294489_wp,      -2.853846527970828E-002_wp, -4.716562899050264E-003_wp, &
    -8.272119708708751E-004_wp, -1.490290753800987E-004_wp, -2.723995785749091E-005_wp, &
    -5.022464331565194E-006_wp, -9.312681555797506E-007_wp, -1.733446482843680E-007_wp, &
    -3.235539170874715E-008_wp &
    ]

    frc_tab = [&
    -0.136218726008017_wp,      -9.523372817110675E-003_wp, -9.972161552254323E-004_wp, &
    -1.283769347101531E-004_wp, -1.829884750086477E-005_wp, -2.769395715511576E-006_wp, &
    -4.358969043562504E-007_wp, -7.052837906523818E-008_wp, -1.164694761021766E-008_wp, &
    -1.953803439326774E-009_wp &
    ]

    electro%key = ELECTROSTATIC_SPME
    electro%smear = SMEARING_SLATER_TRUNCATED
    electro%r_smear = 2.0_wp
    electro%b_smear = 0

    Call electro%init_erf_tables(10)
    Call electro%erfcgen(10.0_wp, 1.0_wp)

    Call assert(electro%erfc_over_r%table, pot_tab, "Slater (apprx) smearing tabulated energies differ from expected", &
                passed_accum=passed)
    Call assert(electro%erfc_gamma%table, frc_tab, "Slater (apprx) smearing tabulated forces differ from expected", &
                passed_accum=passed)

    smear = slater_apprx_smearing(1.0_wp, 2.0_wp, 0)
    Call assert(smear%energy, 0.551819161757164_wp, &
              "Slater (apprx) smearing direct energy differs from expected", passed_accum=passed)
    Call assert(smear%force, 0.9196986029286066_wp, &
              "Slater (apprx) smearing direct force differs from expected", passed_accum=passed)

  End Subroutine slater_apprx_smearing_test

  Subroutine gaussian_smearing_test(passed)
    Logical, Intent(InOut) :: passed
    Type(electrostatic_type) :: electro
    Real(Kind=wp), Dimension(10) :: pot_tab, frc_tab

    pot_tab = [&
    -0.322360565544624_wp,      -7.157715649568776E-002_wp, -1.541995601110098E-002_wp, &
    -2.763313956707276E-003_wp, -3.859639528612339E-004_wp, -4.070354633940066E-005_wp, &
    -3.182422461616815E-006_wp, -1.823495234085744E-007_wp, -7.597671202102616E-009_wp, &
    -2.289177619965212E-010_wp &
    ]
    frc_tab = [&
    -0.176161482357041_wp,      -1.911823481824899E-002_wp, -2.982009464290113E-003_wp, &
    -4.568172199638683E-004_wp, -5.849970142689654E-005_wp, -5.852746039275783E-006_wp, &
    -4.421577131798989E-007_wp, -2.474072821432637E-008_wp, -1.013136528933572E-009_wp, &
    -3.012377964642021E-011_wp &
    ]

    electro%key = ELECTROSTATIC_SPME
    electro%smear = SMEARING_GAUSSIAN
    electro%r_smear = 2.0_wp
    Call electro%init_erf_tables(10)
    Call electro%erfcgen(10.0_wp, 1.0_wp)

    Call assert(electro%erfc_over_r%table, pot_tab, "Gaussian smearing tabulated energies differ from expected", &
                passed_accum=passed)
    Call assert(electro%erfc_gamma%table, frc_tab, "Gaussian smearing tabulated forces differ from expected", &
                passed_accum=passed)

  End Subroutine gaussian_smearing_test

End Module
