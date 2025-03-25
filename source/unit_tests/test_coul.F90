Module test_coul 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_5 module for unit testing the point charge  Coulombic 
  ! interactions available in DL_POLY_5 including possible corrections 
  ! for charge smearing.
  !
  ! copyright - daresbury laboratory
  ! author    - b.t.speake march 2025
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use asserts,            Only: assert 
  Use configuration,      Only: configuration_type
  Use errors_warnings,    Only: error,error_alloc,error_dealloc
  Use kinds,              Only: wp, wi
  Use electrostatic,      Only: electrostatic_type, ELECTROSTATIC_SPME, SMEARING_LINEAR, & 
                                SMEARING_SLATER_EXP, SMEARING_SLATER_TRUNCATED, &
                                SMEARING_GAUSSIAN, ELECTROSTATIC_COULOMB, ELECTROSTATIC_DDDP, &
                                ELECTROSTATIC_COULOMB_FORCE_SHIFT, ELECTROSTATIC_COULOMB_REACTION_FIELD
  Use neighbours,         Only: neighbours_type
  Use statistics,         Only: stats_type
  Use coul_spole,         Only: coul_cp_forces, coul_dddp_forces, coul_fscp_forces, coul_rfp_forces

Implicit None 

Contains 

  Subroutine run_coulomb_tests(passed)
    Logical, Intent(Out) :: passed

    passed = .true. 
    Call test_direct_coulomb(passed)
    Call test_DDDP(passed)
    Call test_force_shifted(passed)
    Call test_force_shifted_damped(passed)
    Call test_reaction_field(passed)
    Call test_reaction_field_damped(passed)
  End Subroutine run_coulomb_tests

  Subroutine test_direct_coulomb(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [27787.0967000000_wp,  27787.0967000000_wp,  27787.0808938993_wp, &     
               27787.0919834371_wp,  27787.0812800440_wp]
    frc_tab = [-5557.41934000000_wp, -5557.41934000000_wp, -5557.40593148547_wp, &     
               -5557.41435391922_wp, -5557.41862122002_wp]

    electro%key = ELECTROSTATIC_COULOMB 
    electro%r_smear =  2.0_wp 

    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call coul_cp_forces(1, electro, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    Call coul_cp_forces(1, electro, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    Call coul_cp_forces(1, electro, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    Call coul_cp_forces(1, electro, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    Call coul_cp_forces(1, electro, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, config)
    frc(5) = config%parts(1)%fxx    

    Call assert(eng, eng_tab, "Direct coulomb energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "Direct coulomb forces differ from expected values", passed_accum=passed)

  End Subroutine test_direct_coulomb

  Subroutine test_DDDP(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [ 5557.41934000000_wp,  5557.41934000000_wp,  5557.41617877987_wp, &     
                5557.41839668742_wp,  5557.41625600880_wp  ]
    frc_tab = [-2222.96773600000_wp, -2222.96773600000_wp, -2222.96505429709_wp, &     
               -2222.96673878384_wp, -2222.96759224400_wp  ]

    electro%key = ELECTROSTATIC_DDDP
    electro%r_smear =  2.0_wp 

    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call coul_dddp_forces(1, electro, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    Call coul_dddp_forces(1, electro, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    Call coul_dddp_forces(1, electro, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    Call coul_dddp_forces(1, electro, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    Call coul_dddp_forces(1, electro, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, config)
    frc(5) = config%parts(1)%fxx     

    Call assert(eng, eng_tab, "DDDP energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "DDDP forces differ from expected values", passed_accum=passed)

  End Subroutine test_DDDP

  Subroutine test_force_shifted(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [6946.77417500000_wp,  6946.77417500000_wp,  4851.74736118974_wp, &     
               6306.92999804070_wp,  18738.6011860444_wp ]
    frc_tab = [-4168.06450500000_wp, -4168.06450500000_wp, -2331.33336272771_wp, &     
               -3479.16863448204_wp, -4882.24004857187_wp]

    electro%key = ELECTROSTATIC_COULOMB_FORCE_SHIFT
    electro%r_smear =  2.0_wp 

    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, electro, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, electro, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, electro, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, electro, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, electro, config)
    frc(5) = config%parts(1)%fxx  

    Call assert(eng, eng_tab, "Force shifted coulomb energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "Force shifted coulomb forces differ from expected values", passed_accum=passed)

  End Subroutine test_force_shifted

  Subroutine test_force_shifted_damped(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [ 4.300934185702100E-008_wp,  4.300934185702100E-008_wp,  -2035.15659628256_wp, &     
               -632.274983485388_wp,       -2096.06617202174_wp  ]
    frc_tab = [-4.440497520415766E-007_wp, -4.440497520415766E-007_wp,   1836.73114182824_wp, &     
                688.895870073908_wp,        2063.40309260592_wp   ]

    electro%key = ELECTROSTATIC_COULOMB_FORCE_SHIFT
    electro%r_smear =  2.0_wp 
    electro%damp = .true.
    electro%damping = 1.0_wp 
    
    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, electro, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, electro, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, electro, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, electro, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, electro, config)
    frc(5) = config%parts(1)%fxx  

    Call assert(eng, eng_tab, "Force shifted coulomb (w/ fennel damping) energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "Force shifted coulomb (w/ fennel damping) forces differ from expected values", passed_accum=passed)

  End Subroutine test_force_shifted_damped

  Subroutine test_reaction_field(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [18524.7311333333_wp,  18524.7311333333_wp, 17060.7123075641_wp, &     
               18087.8658354129_wp, 17096.4784373373_wp ]
    frc_tab = [-3704.94622666667_wp,  -3704.94622666667_wp, -2463.00059400702_wp, &     
               -3243.11719743645_wp,  -2323.92313729489_wp]

    electro%key = ELECTROSTATIC_COULOMB_REACTION_FIELD
    electro%r_smear =  2.0_wp
    electro%eps = 1.5_wp

    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, electro, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, electro, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, electro, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, electro, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    electro%initialised = .false.
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, electro, config)
    frc(5) = config%parts(1)%fxx  

    Call assert(eng, eng_tab, "Reaction field coulomb energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "Reaction field coulomb forces differ from expected values", passed_accum=passed)

  End Subroutine test_reaction_field

  Subroutine test_reaction_field_damped(passed)
    Logical, Intent(InOut) :: passed

    Type(electrostatic_type)                 :: electro
    Type(neighbours_type)                    :: neigh
    Type(stats_type)                         :: stats
    Type(configuration_type)                 :: config

    Real(Kind=wp) :: vir 
    Real(kind=wp), Allocatable, Dimension(:)    :: rrt, xxt, yyt, zzt
    Integer(Kind=wi) :: j, k, fail
    Real(Kind=wp), Dimension(5) :: eng = 0.0_wp 
    Real(Kind=wp), Dimension(5) :: frc = 0.0_wp
    Real(Kind=wp), Dimension(5) :: eng_tab, frc_tab

    ! energy and force tables (pure, linear, slater, slater approx, gaussian)
    eng_tab = [   2.867289457134733E-008_wp, 2.867289457134733E-008_wp, -1356.77106418838_wp, &     
               -421.516655656926_wp, -1397.37744801450_wp   ]
    frc_tab = [ -2.960331680277178E-007_wp,  -2.960331680277178E-007_wp,   1224.48742788549_wp, &     
               459.263913382605_wp, 1375.60206173728_wp]

    electro%key = ELECTROSTATIC_COULOMB_REACTION_FIELD
    electro%r_smear =  2.0_wp 
    electro%eps = 1.5_wp 
    electro%damp = .true.
    electro%damping = 1.0_wp 
    
    Call setup_fake_system(neigh, stats, config)

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    ! calculate interatomic distances
    Do k = 1, neigh%list(0, 1)
      j = neigh%list(k, 1)
      xxt(k) = config%parts(1)%xxx - config%parts(j)%xxx
      yyt(k) = config%parts(1)%yyy - config%parts(j)%yyy
      zzt(k) = config%parts(1)%zzz - config%parts(j)%zzz
      rrt(k) = Sqrt(xxt(k)*xxt(k)+yyt(k)*yyt(k)+zzt(k)*zzt(k))
    End Do

    ! (1) pure Coulomb interactions
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(1), vir, stats, neigh, electro, config)
    frc(1) = config%parts(1)%fxx 
    
    ! (2) linear charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_LINEAR 
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(2), vir, stats, neigh, electro, config)
    frc(2) = config%parts(1)%fxx 

    ! (3) slater (exp) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_EXP  
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(3), vir, stats, neigh, electro, config)
    frc(3) = config%parts(1)%fxx 

    ! (4) slater (truncated) charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_SLATER_TRUNCATED 
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(4), vir, stats, neigh, electro, config)
    frc(4) = config%parts(1)%fxx 

    ! (5) Gaussian charge smearing 
    config%parts(:)%fxx = 0.0_wp
    electro%smear = SMEARING_GAUSSIAN
    Call init_electro_tables(electro, neigh%cutoff)
    Call coul_fscp_forces(1, xxt, yyt, zzt, rrt, eng(5), vir, stats, neigh, electro, config)
    frc(5) = config%parts(1)%fxx  

    Call assert(eng, eng_tab, "Force shifted coulomb (w/ fennel damping) energies differ from expected values", passed_accum=passed)
    Call assert(frc, frc_tab, "Force shifted coulomb (w/ fennel damping) forces differ from expected values", passed_accum=passed)

  End Subroutine test_reaction_field_damped

  Subroutine setup_fake_system(neigh, stats, config)
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(stats_type),         Intent(InOut) :: stats
    Type(configuration_type), Intent(InOut) :: config

    neigh%max_list = 2
    neigh%max_exclude = 0
    neigh%cutoff = 10.0_wp
    call neigh%init_list(2)
    neigh%list(0, 1) = 1
    neigh%list(1, 1) = 2

    config%mxatms = 2
    config%mxatdm = 2
    config%natms = 10
    Call config%init_read()
    Call config%init()
    config%parts(2)%xxx = 5.0_wp
    config%parts(1)%chge = 1.0_wp 
    config%parts(2)%chge = 1.0_wp
    config%ltype = 1

    stats%collect_pp_eng_str = .false.
    stats%stress = 0.0_wp
    config%parts(:)%fxx = 0.0_wp
    config%parts(:)%fyy = 0.0_wp
    config%parts(:)%fzz = 0.0_wp


  End Subroutine setup_fake_system

  Subroutine init_electro_tables(electro, rcut)

    Type(electrostatic_type), Intent(InOut) :: electro 
    Real(Kind=wp),            Intent(In   ) :: rcut 
    Integer(Kind=wi)                        :: fail
 
    electro%initialised = .false.
    If (electro%erfc%initialised) Then 
      Deallocate (electro%erfc%table, stat=fail)
      If (fail > 0) Call error_dealloc('electro%erfc%table', 'init_interp_table')
      electro%erfc%spacing = -1.0_wp
      electro%erfc%recip_spacing = -1.0_wp
      electro%erfc%nsamples = 0
      electro%erfc%initialised = .false.
    End If 
    If (electro%erfc_deriv%initialised) Then 
      Deallocate (electro%erfc_deriv%table, stat=fail)
      If (fail > 0) Call error_dealloc('electro%erfc_deriv%table', 'init_interp_table')
      electro%erfc_deriv%spacing = -1.0_wp
      electro%erfc_deriv%recip_spacing = -1.0_wp
      electro%erfc_deriv%nsamples = 0
      electro%erfc_deriv%initialised = .false.
    End If 
    Call electro%init_erf_tables(Max(1004, Nint(rcut / 0.01_wp) + 4))
    Call electro%erfcgen(rcut, electro%damping)

  End Subroutine init_electro_tables

End Module test_coul 
