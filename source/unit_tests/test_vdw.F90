Module test_vdw

  Use asserts,                  Only: assert
  Use constants,                Only: wp
  Use vdw,                      Only: vdw_type, NUM_VDW_POTS, vdw_forces_direct
  Use neighbours,               Only: neighbours_type
  Use statistics,               Only: stats_type
  Use configuration,            Only: configuration_type
  Use two_body_potentials,      Only: potential, &
                                      potential_energy, &
                                      potential_holder, &
                                      LJ, lj_coh, LJ126, n_m, &
                                      nm_shift, morse, morse12, &
                                      buckingham, bhm, hbond, &
                                      wca, dpd, ndpd, amoeba, &
                                      rydberg, zbl, zblb, fm, zbls, &
                                      sanderson, MDF, ljf, mlj, &
                                      mbuck, mlj126, sw

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
        4.0000000000000000_wp, -0.8948393168143698_wp, -0.4921875000000000_wp, -1765.997274794514_wp]
    Real(kind=wp), Dimension(NUM_VDW_POTS), Parameter :: expected_v = [&
        0.0000000000000000_wp, -195072.00000000000_wp, -288.00000000000000_wp, 17.696734670143684_wp, &
        -45.196300066288472_wp, 8.0000000000000000_wp, -16899.173553719396_wp, -2300.0595394172847_wp, &
        12.000000000000000_wp, -0.5000000000000000_wp, -3719.0014773362891_wp, -192000.00000000000_wp, &
        -2348.0595394172847_wp, 0.71653131057378916_wp, -50044.066263312896_wp, -5.2445617006037446E+031_wp, &
        -36135.103419546060_wp, -195072.00000000000_wp, 17.696734670143684_wp, 0.0000000000000000_wp, &
        -40.000000000000000_wp, -0.1988531815143044_wp, -0.09375000000000000_wp, 15903.665444480088_wp]


    Type(vdw_type) :: test
    Type(neighbours_type) :: neigh
    Type(stats_type) :: stats
    Type(configuration_type) :: config
    Real(wp), Dimension(1:NUM_VDW_POTS) :: eng = 0.0_wp, gamma = 0.0_wp
    Class(potential_holder), Allocatable :: pots(:)
    Integer :: i

    passed = .true.

    Call setup_fake_system(test, pots, neigh, stats, config)

    do i = 1, NUM_VDW_POTS
      ! Cycle through each pot
      test%ltp = i
      test%potentials(1) = pots(i)
      Call test%update_potential_cutoffs()
      
      stats%stress = 0.0_wp
      config%parts(:)%fxx = 0.0_wp
      config%parts(:)%fyy = 0.0_wp
      config%parts(:)%fzz = 0.0_wp

      Call vdw_forces_direct(1, ones, ones, ones, ones, eng(i), gamma(i), stats, neigh, test, config)
    end do

    Call assert(eng, expected_e, "VdW Energies differ from expected", passed_accum = passed)
    Call assert(gamma, expected_v, "VdW Virials differ from expected", passed_accum = passed)

  end Subroutine test_forces_direct

  Subroutine setup_fake_system(test, pots, neigh, stats, config)
    Type(vdw_type), Intent(InOut) :: test
    Class(potential_holder), Allocatable, Intent(InOut) :: pots(:)
    Type(neighbours_type), Intent(InOut) :: neigh
    Type(stats_type), Intent(InOut) :: stats
    Type(configuration_type), Intent(InOut) :: config

    test%max_vdw = 1
    test%max_param = 7
    test%newjob = .false.
    test%cutoff = 10.0_wp
    Call test%init()
    Call test%init_direct()
    test%param(:,1) = params
    test%list = 1

    Allocate(pots(1:NUM_VDW_POTS+1))

    Allocate(LJ126::pots(1)%p)
    Call pots(1)%p%set_parameters(params)
    Allocate(LJ::pots(2)%p)
    Call pots(2)%p%set_parameters(params)
    Allocate(n_m::pots(3)%p)
    Call pots(3)%p%set_parameters(params)
    Allocate(buckingham::pots(4)%p)
    Call pots(4)%p%set_parameters(params)
    Allocate(bhm::pots(5)%p)
    Call pots(5)%p%set_parameters(params)
    Allocate(hbond::pots(6)%p)
    Call pots(6)%p%set_parameters(params)
    Allocate(nm_shift::pots(7)%p)
    Call pots(7)%p%set_parameters(params)
    Allocate(morse::pots(8)%p)
    Call pots(8)%p%set_parameters(params)
    Allocate(wca::pots(9)%p)
    Call pots(9)%p%set_parameters(params)
    Allocate(dpd::pots(10)%p)
    Call pots(10)%p%set_parameters(params)
    Allocate(amoeba::pots(11)%p)
    Call pots(11)%p%set_parameters(params)
    Allocate(lj_coh::pots(12)%p)
    Call pots(12)%p%set_parameters(params)
    Allocate(morse12::pots(13)%p)
    Call pots(13)%p%set_parameters(params)
    Allocate(rydberg::pots(14)%p)
    Call pots(14)%p%set_parameters(params)
    Allocate(zbl::pots(15)%p)
    Call pots(15)%p%set_parameters(params)
    Allocate(zbls::pots(16)%p)
    Call pots(16)%p%set_parameters(params)
    Allocate(zblb::pots(17)%p)
    Call pots(17)%p%set_parameters(params)
    Allocate(mlj::pots(18)%p)
    Call pots(18)%p%set_parameters(params)
    Allocate(mbuck::pots(19)%p)
    Call pots(19)%p%set_parameters(params)
    Allocate(mlj126::pots(20)%p)
    Call pots(20)%p%set_parameters(params)
    Allocate(ljf::pots(21)%p)
    Call pots(21)%p%set_parameters(params)
    Allocate(sanderson::pots(22)%p)
    Call pots(22)%p%set_parameters(params)
    Allocate(ndpd::pots(23)%p)
    Call pots(23)%p%set_parameters(params)
    Allocate(sw::pots(24)%p)
    Call pots(24)%p%set_parameters(params)

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
