Module test_vdw

  Use asserts, only : assert
  Use constants, only : wp
  Use vdw, only : vdw_type, NUM_VDW_POTS, vdw_forces_direct! , LJ, lj_coh, LJ126, &
       ! n_m, nm_shift, morse, morse_12, buckingham, bhm, hbond, &
       ! wca, dpd, amoeba, rydberg, zbl, fm, zbls, zblb, MDF, ljf, &
       ! mlj, mbuck, mlj126
  Use neighbours, only : neighbours_type
  Use statistics, only : stats_type
  Use configuration, only : configuration_type

  Implicit None

  Private

  Real(kind=wp), Dimension(2), Parameter :: ones = 1.0_wp
  Real(kind=wp), Parameter :: r = 1.0_wp
  Integer :: tmp
  Real(kind=wp), Dimension(7), Parameter :: params = [(real(tmp,wp), tmp=1,7)]

  Public :: run_vdw_tests

Contains

  Subroutine run_vdw_tests(passed)
      Logical, Intent(  Out) :: passed

      passed = .true.
    ! Call test_indiv_funcs()
    Call test_forces_direct(passed)

  end Subroutine run_vdw_tests

  Subroutine test_forces_direct(passed)

    Logical, Intent(  Out) :: passed

    Real(kind=wp), Dimension(NUM_VDW_POTS), Parameter :: expected_e =  [&
        -1.0000000000000000_wp, 16128.000000000000_wp, 80.000000000000000_wp, -2.3934693402873668_wp, &
        45.598150033144236_wp, -1.0000000000000000_wp, 4567.8016528926664_wp, 363.25771964635982_wp, &
        1.0000000000000000_wp, 0.25000000000000000_wp, 555.50928442334305_wp, 15616.000000000000_wp, &
        367.25771964635982_wp, 2.1495939317213679_wp, 16882.652704957309_wp, 3.8142266913481757E+030_wp, &
        11761.890149027380_wp, 16128.000000000000_wp, -2.3934693402873668_wp, -1.0000000000000000_wp, &
        4.0000000000000000_wp]
    Real(kind=wp), Dimension(NUM_VDW_POTS), Parameter :: expected_v = [&
        0.0000000000000000_wp, -195072.00000000000_wp, -288.00000000000000_wp, 17.696734670143684_wp, &
        -45.196300066288472_wp, 8.0000000000000000_wp, -16899.173553719396_wp, -2300.0595394172847_wp, &
        12.000000000000000_wp, -0.5000000000000000_wp, -3719.0014773362891_wp, -192000.00000000000_wp, &
        -2348.0595394172847_wp, 0.71653131057378916_wp, -50044.066263312896_wp, -5.2445617006037446E+031_wp, &
        -36135.103419546060_wp, -195072.00000000000_wp, 17.696734670143684_wp, 0.0000000000000000_wp, &
        -40.000000000000000_wp]


    Type(vdw_type) :: test
    Type(neighbours_type) :: neigh
    Type(stats_type) :: stats
    Type(configuration_type) :: config
    Real(wp), Dimension(NUM_VDW_POTS) :: eng = 0.0_wp, gamma = 0.0_wp
    Integer :: i

    passed = .true.

    Call setup_fake_system(neigh, stats, config)

    test%max_vdw = 1
    test%max_param = 7
    test%newjob = .false.
    test%cutoff = 10.0_wp
    Call test%init()
    Call test%init_direct()
    test%param(:,1) = params
    test%list = 1

    do i = 1, NUM_VDW_POTS
      ! Cycle through each pot
      test%ltp = i
      stats%stress = 0.0_wp
      config%parts(:)%fxx = 0.0_wp
      config%parts(:)%fyy = 0.0_wp
      config%parts(:)%fzz = 0.0_wp

      Call vdw_forces_direct(1, ones, ones, ones, ones, eng(i), gamma(i), stats, neigh, test, config)
    end do

    Call assert(eng, expected_e, "VdW Energies differ from expected", passed_accum = passed)
    Call assert(gamma, expected_v, "VdW Virials differ from expected", passed_accum = passed)

  end Subroutine test_forces_direct

  ! Subroutine test_indiv_funcs()

  !   Real(kind=wp), Dimension(NUM_VDW_POTS), Parameter :: expected_e =  [&
  !        16128.000000000000, 15616.000000000000, -1.0000000000000000, 80.000000000000000, &
  !        4567.8016528925955, 363.25771964635982, 367.25771964635982, -2.3934693402873668, &
  !        45.598150033144236, -1.0000000000000000, 1.0000000000000000, 0.12500000000000000, &
  !        555.50928442334305, 2.1495939317213679, 16882.652704957312, 3.8142266913481757E+030, &
  !        11761.890149027382, 4.0000000000000000, 16128.000000000000, -2.3934693402873668, &
  !        -1.0000000000000000]
  !   Real(kind=wp), Dimension(NUM_VDW_POTS), Parameter :: expected_v = [&
  !        195072.00000000000, 192000.00000000000, 0.0000000000000000, 288.00000000000000, &
  !        16899.173553719134, 2300.0595394172847, 2348.0595394172847, -17.696734670143684, &
  !        45.196300066288472, -8.0000000000000000, -12.000000000000000, 0.12500000000000000, &
  !        3719.0014773362891, -0.71653131057378916, 50044.066263312910, 5.2445617006037455E+031, &
  !        36135.103419546067, 40.000000000000000, 195072.00000000000, -17.696734670143684, &
  !        0.0000000000000000]

  !   Real(kind=wp), Dimension(NUM_VDW_POTS) :: eng, gamma

  !   Call LJ         (r, params, eng( 1 ), gamma( 1 ))
  !   Call lj_coh     (r, params, eng( 2 ), gamma( 2 ))
  !   Call LJ126      (r, params, eng( 3 ), gamma( 3 ))
  !   Call n_m        (r, params, eng( 4 ), gamma( 4 ))
  !   Call nm_shift   (r, params, eng( 5 ), gamma( 5 ))
  !   Call morse      (r, params, eng( 6 ), gamma( 6 ))
  !   Call morse_12   (r, params, eng( 7 ), gamma( 7 ))
  !   Call buckingham (r, params, eng( 8 ), gamma( 8 ))
  !   Call bhm        (r, params, eng( 9 ), gamma( 9 ))
  !   Call hbond      (r, params, eng(10 ), gamma(10 ))
  !   Call wca        (r, params, eng(11 ), gamma(11 ))
  !   Call dpd        (r, params, eng(12 ), gamma(12 ))
  !   Call amoeba     (r, params, eng(13 ), gamma(13 ))
  !   Call rydberg    (r, params, eng(14 ), gamma(14 ))
  !   Call zbl        (r, params, eng(15 ), gamma(15 ))
  !   Call zbls       (r, params, eng(16 ), gamma(16 ))
  !   Call zblb       (r, params, eng(17 ), gamma(17 ))
  !   Call ljf        (r, params, eng(18 ), gamma(18 ))
  !   Call mlj        (r, params, eng(19 ), gamma(19 ))
  !   Call mbuck      (r, params, eng(20 ), gamma(20 ))
  !   Call mlj126     (r, params, eng(21 ), gamma(21 ))

  !   Call assert(all(abs(eng - expected_e) < 1.0e-6), "VdW Energies differ from expected")
  !   Call assert(all(abs(gamma - expected_v) < 1.0e-6), "VdW Virials differ from expected")

  ! end Subroutine test_indiv_funcs

  Subroutine setup_fake_system(neigh, stats, config)
    Type(neighbours_type), Intent(InOut) :: neigh
    Type(stats_type), Intent(InOut) :: stats
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
    config%parts(2)%xxx = 1.0_wp
    config%ltype = 1

    stats%collect_pp = .false.

  end Subroutine setup_fake_system


end Module test_vdw
