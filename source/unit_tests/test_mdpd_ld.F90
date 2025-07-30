Module test_mdpd_ld 

  Use vdw,                 Only: vdw_type, VDW_MDPD
  Use configuration,       Only: configuration_type
  Use domains,             Only: domains_type
  Use neighbours,          Only: neighbours_type
  Use comms,               Only: comms_type
  Use asserts,             Only: assert
  Use constants,           Only: wp 
  Use two_body_potentials, Only: potential_holder

  Implicit None
  
  Real(kind=wp), Dimension(2), Parameter :: ones = 1.0_wp
  Real(kind=wp),               Parameter :: r = 1.0_wp
  Integer                                :: tmp
  Real(kind=wp), Dimension(7), Parameter :: params = [(real(tmp,wp), tmp=1,7)]

  Public :: run_mdpd_tests

Contains 

  Subroutine run_mdpd_tests(passed) 
    Logical, Intent(Out) :: passed
    passed = .true.
    Call mdpd_ld_test(passed)
  End Subroutine run_mdpd_tests

  Subroutine mdpd_ld_test(passed)
    Use mdpd, Only: weight_rho
    Logical, Intent(InOut) :: passed

    Type(vdw_type)                       :: vdws 
    Type(configuration_type)             :: config
    Type(neighbours_type)                :: neigh
    Class(potential_holder), Allocatable :: pots(:)
    Integer                              :: i, j, k, aj, ai, key
    Real(Kind=wp)                        :: rho, xxt, yyt, zzt, rrr

    Call setup_fake_system(vdws, pots, neigh, config)

    vdws%ltp = VDW_MDPD 
    vdws%potentials(1) = pots(1) 
    vdws%l_direct = .true.
    Call vdws%set_constants() 

    config%parts(:)%fxx = 0.0_wp 
    config%parts(:)%fyy = 0.0_wp 
    config%parts(:)%fzz = 0.0_wp 

    Allocate(vdws%mdpd_params%rho(config%mxatms))
    vdws%mdpd_params%rho = 0.0_wp 

    Do i = 1, config%natms 
      Do k = 1, neigh%list(0, i) 
        j = neigh%list(k, i) 
        xxt = config%parts(i)%xxx - config%parts(j)%xxx
        yyt = config%parts(i)%yyy - config%parts(j)%yyy
        zzt = config%parts(i)%zzz - config%parts(j)%zzz
        rrr = Sqrt(xxt**2 + yyt**2 + zzt**2)

        ai = config%ltype(i)
        aj = config%ltype(j)

        If (ai > aj) Then 
          key = ai * (ai - 1) / 2 + aj
        Else 
          key = aj * (aj - 1) / 2 + ai
        End If

        rho = weight_rho(rrr, vdws%mdpd_params%rc(key), vdws%mdpd_params%rd(key), vdws%mdpd_params%n(key))
        vdws%mdpd_params%rho(i) = vdws%mdpd_params%rho(i) + rho  
        vdws%mdpd_params%rho(j) = vdws%mdpd_params%rho(j) + rho
      End Do 
    End Do 

    Call assert(vdws%mdpd_params%rho(1), 0.0019337325585665284_wp, passed_accum = passed) 
    Call assert(vdws%mdpd_params%rho(2), 0.0019337325585665284_wp, passed_accum = passed)

  End Subroutine 

  Subroutine setup_fake_system(test, pots, neigh, config)
    Use two_body_potentials, Only: mdpd

    Type(vdw_type),                       Intent(InOut) :: test
    Class(potential_holder), Allocatable, Intent(InOut) :: pots(:)
    Type(neighbours_type),                Intent(InOut) :: neigh
    Type(configuration_type),             Intent(InOut) :: config

    test%max_vdw = 1
    test%max_param = 7
    test%cutoff = 10.0_wp
    Call test%init()
    Call test%init_direct()
    test%param(:,1) = params
    test%list = 1

    Allocate(pots(1))
    Allocate(mdpd::pots(1)%p)
    Call pots(1)%p%set_parameters(params)
    Call setup_mdpd_parameters(test)

    neigh%max_list = 2
    neigh%max_exclude = 0
    neigh%cutoff = 10.0_wp
    call neigh%init_list(2)
    neigh%list(0, 1) = 1
    neigh%list(1, 1) = 2

    config%mxatms = 2
    config%mxatdm = 2
    config%natms = 2
    Call config%init_read()
    Call config%init()
    config%parts(2)%xxx = 1.0_wp
    config%ltype = 1

  end Subroutine setup_fake_system

  Subroutine setup_mdpd_parameters(vdws) 
    Type(vdw_type), Intent(InOut) :: vdws 
    vdws%mdpd_params%b = 1 
    
    Allocate(vdws%mdpd_params%rd(1))
    Allocate(vdws%mdpd_params%m(1))
    Allocate(vdws%mdpd_params%n(1))
    Allocate(vdws%mdpd_params%rc(1))
    ! Allocate(vdws%mdpd_params%rho(2))

    vdws%mdpd_params%n = 2.0_wp 
    vdws%mdpd_params%m = 2.0_wp 
    vdws%mdpd_params%rd = 10.0_wp 
    vdws%mdpd_params%rc = 1.0_wp
    ! vdws%mdpd_params%rho = 1.0_wp 
  End Subroutine 

End Module test_mdpd_ld
