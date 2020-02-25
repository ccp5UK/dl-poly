Module nvt_langevin
  Use comms,           Only: comms_type
  Use configuration,   Only: configuration_type
  Use constants,       Only: zero_plus
  Use constraints,     Only: apply_rattle,&
                             apply_shake,&
                             constraints_tags,&
                             constraints_type
  Use core_shell,      Only: core_shell_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info
  Use kinds,           Only: wp
  Use kinetics,        Only: getknr,&
                             getvom,&
                             kinstresf,&
                             kinstress,&
                             kinstrest
  Use langevin,        Only: langevin_forces
  Use numerics,        Only: images,&
                             seed_type
  Use pmf,             Only: pmf_tags,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             no_squish,&
                             rigid_bodies_stress,&
                             rigid_bodies_type
  Use shared_units,    Only: update_shared_units
  Use statistics,      Only: stats_type
  Use thermostat,      Only: VV_FIRST_STAGE,&
                             adjust_timestep,&
                             thermostat_type
  Use timer,           Only: timer_type
  Use ttm,             Only: eltemp_max,eltemp_min,&
                             ttm_type
  Use ttm_utils,       Only: Gep,&
                             calcchies

  Implicit None

  Private

  Public :: nvt_l0_vv, nvt_l1_vv, nvt_l2_vv

Contains

  Subroutine nvt_l0_vv(stage, lvar, mndis, mxdis, mxstp, tstep, nstep, strkin, engke, &
                       cshell, cons, pmf, stat, thermo, domain, tmr, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with Langevin thermostat, based
    ! on Langevin impuls (LI) integration
    !
    ! Ref: Jesus A. Izaguirre, `Langevin Stabilisation of Multiscale
    !      Mollified Dynamics', editors: A. Brandt, K. Binder, J. Bernholc,
    !      Multiscale Computational Methods in Chemistry and Physics,
    !      January 2001, vol. 117 of NATO Science Series: Series III -
    !      Computer and System Sciences, pp. 34-47, IOS Press, Amsterdam
    !
    ! (brownian dynamics is not symplectic due to the random forces)
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), engke
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:9), i
    Logical                    :: safe
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: hstep, rstep, scl, scl1, scr, scr1, scv, scv1, t0, t1, t2, tmp, &
                                  vom(1:3)
    Real(Kind=wp), Allocatable :: fxl(:), fxr(:), fxt(:), fyl(:), fyr(:), fyt(:), fzl(:), fzr(:), &
                                  fzt(:), oxt(:), oyt(:), ozt(:), vxt(:), vyt(:), vzt(:), xxt(:), &
                                  yyt(:), zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms), oyt(1:config%mxatms), ozt(1:config%mxatms), Stat=fail(6))
    End If
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(7))
    Allocate (vxt(1:config%mxatms), vyt(1:config%mxatms), vzt(1:config%mxatms), Stat=fail(8))
    Allocate (fxt(1:config%mxatms), fyt(1:config%mxatms), fzt(1:config%mxatms), Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l0 allocation failure'
      Call error(0, message)
    End If

    If (thermo%newjob_0) Then
      thermo%newjob_0 = .false.

      ! set number of constraint+pmf shake iterations
      If (cons%megcon > 0 .or. pmf%megpmf > 0) thermo%mxkit = 1
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms
      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms
      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
      End If
    End If

    ! timestep derivatives
    hstep = 0.5_wp * tstep
    rstep = 1.0_wp / tstep

    ! first pass of velocity verlet algorithm
    If (stage == VV_FIRST_STAGE) Then

      ! store initial values
      Do i = 1, config%natms
        xxt(i) = config%parts(i)%xxx
        yyt(i) = config%parts(i)%yyy
        zzt(i) = config%parts(i)%zzz

        vxt(i) = config%vxx(i)
        vyt(i) = config%vyy(i)
        vzt(i) = config%vzz(i)

        fxt(i) = config%parts(i)%fxx
        fyt(i) = config%parts(i)%fyy
        fzt(i) = config%parts(i)%fzz
      End Do

      ! Set afresh Langevin random forces

      Allocate (fxr(1:config%mxatms), fyr(1:config%mxatms), fzr(1:config%mxatms), Stat=fail(1))
      Allocate (fxl(1:config%mxatms), fyl(1:config%mxatms), fzl(1:config%mxatms), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'nvt_l0 allocation failure+'
        Call error(0, message)
      End If

      Call langevin_forces(nstep, thermo%temp, tstep, thermo%chi, fxr, fyr, fzr, cshell, config, seed)
      Call langevin_forces(-nstep, thermo%temp, tstep, thermo%chi, fxl, fyl, fzl, cshell, config, seed)

      100 Continue

      ! constraint virial and stress tensor

      If (cons%megcon > 0) Then
        stat%vircon = 0.0_wp
        stat%strcon = 0.0_wp
      End If

      ! PMF virial and stress tensor

      If (pmf%megpmf > 0) Then
        stat%virpmf = 0.0_wp
        stat%strpmf = 0.0_wp
      End If

      ! Create primitive scalers and adjust/increase timestep if need be
      ! when Cholesky factorisation is compromised

      t0 = Exp(-thermo%chi * tstep)
      t1 = (1.0_wp - t0) / (thermo%chi)
      t2 = (1.0_wp - t0**2) / (2 * thermo%chi)

      safe = .true.
      Do
        tmp = t1**2 / t2
        If (tstep - tmp >= zero_plus) Then
          If ((.not. safe)) Then
            Write (message, "('timestep increased due to impossibility of integration, new timestep is:',3x,1p,e16.8,/)") &
              tstep

            Call info(message, .true.)
          Endif
          Exit
        Else
          safe = .false.
          tstep = tmp + 1.0e-10_wp
          t0 = Exp(-thermo%chi * tstep)
          t1 = (1.0_wp - t0) / (thermo%chi)
          t2 = (1.0_wp - t0**2) / (2 * thermo%chi)
        End If
      End Do

      ! Create complex scalers

      scr = (t1 - t2) / Sqrt(t2)
      scl = Sqrt(tstep - (t1**2) / t2)
      scv = Sqrt(t2)

      scr1 = (t1 - t2) / Sqrt(t2 * tstep) / thermo%chi
      scl1 = Sqrt(1.0_wp - (t1**2) / (t2 * tstep)) / thermo%chi
      scv1 = Sqrt(t2 / tstep)

      ! update velocity and position

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then

          ! Half-kick velocity

          tmp = hstep / config%weight(i)
          config%vxx(i) = vxt(i) + tmp * fxt(i)
          config%vyy(i) = vyt(i) + tmp * fyt(i)
          config%vzz(i) = vzt(i) + tmp * fzt(i)

          tmp = tstep / config%weight(i)

          ! Full time fluctuations on positions using half-kick velocity

          config%parts(i)%xxx = xxt(i) + config%vxx(i) * t1 + tmp * (fxr(i) * scr1 + fxl(i) * scl1)
          config%parts(i)%yyy = yyt(i) + config%vyy(i) * t1 + tmp * (fyr(i) * scr1 + fyl(i) * scl1)
          config%parts(i)%zzz = zzt(i) + config%vzz(i) * t1 + tmp * (fzr(i) * scr1 + fzl(i) * scl1)

          ! Full time fluctuations on half-kick velocity

          config%vxx(i) = config%vxx(i) * t0 + tmp * fxr(i) * scv1
          config%vyy(i) = config%vyy(i) * t0 + tmp * fyr(i) * scv1
          config%vzz(i) = config%vzz(i) * t0 + tmp * fzr(i) * scv1

        End If
      End Do

      ! SHAKE procedures

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep, oxt, oyt, ozt, &
                         lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
      End If

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, tmp, comm)) Then
          Call info(message, .true.)

          ! scale Langevin random forces

          Do i = 1, config%natms
            fxr(i) = fxr(i) * tmp
            fyr(i) = fyr(i) * tmp
            fzr(i) = fzr(i) * tmp

            fxl(i) = fxl(i) * tmp
            fyl(i) = fyl(i) * tmp
            fzl(i) = fzl(i) * tmp
          End Do

          ! restart vv1

          Go To 100
        End If
      End If

      Deallocate (fxr, fyr, fzr, Stat=fail(1))
      Deallocate (fxl, fyl, fzl, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'nvt_l0 deallocation failure+'
        Call error(0, message)
      End If

      ! second stage of velocity verlet algorithm

    Else

      ! update velocity (another half-kick)

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
          config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
          config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! remove system centre of mass velocity

      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstress(strkin, config, comm)
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%deallocate_work()
      Deallocate (oxt, oyt, ozt, Stat=fail(6))
    End If
    Deallocate (xxt, yyt, zzt, Stat=fail(7))
    Deallocate (vxt, vyt, vzt, Stat=fail(8))
    Deallocate (fxt, fyt, fzt, Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l0 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nvt_l0_vv

  Subroutine nvt_l1_vv(stage, lvar, mndis, mxdis, mxstp, tstep, nstep, strkin, strknf, &
                       strknt, engke, engrot, strcom, vircom, cshell, cons, pmf, stat, thermo, rigid, &
                       domain, tmr, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet with Langevin thermostat,
    ! based on Langevin impuls (LI) integration
    !
    ! Ref: Jesus A. Izaguirre, `Langevin Stabilisation of Multiscale
    !      Mollified Dynamics', editors: A. Brandt, K. Binder, J. Bernholc,
    !      Multiscale Computational Methods in Chemistry and Physics,
    !      January 2001, vol. 117 of NATO Science Series: Series III -
    !      Computer and System Sciences, pp. 34-47, IOS Press, Amsterdam
    !
    ! (brownian dynamics is not symplectic due to the random forces)
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), strknf(1:9), strknt(1:9), engke, &
                                               engrot, strcom(1:9), vircom
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:14), i, i1, i2, irgd, j, jrgd, krgd, lrgd, matms, rgdtyp
    Logical                    :: safe
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: fmx, fmxl, fmxr, fmy, fmyl, fmyr, fmz, fmzl, fmzr, hstep, mxdr, &
                                  opx, opy, opz, p0, p00, p1, p11, p2, p22, p3, p33, qt0, qt0l, &
                                  qt0r, qt1, qt1l, qt1r, qt2, qt2l, qt2r, qt3, qt3l, qt3r, &
                                  rot(1:9), rstep, scl, scl1, scr, scr1, scv, scv1, t0, t1, t2, &
                                  tmp, tqx, tqxl, tqxr, tqy, tqyl, tqyr, tqz, tqzl, tqzr, trx, &
                                  trxl, trxr, try, tryl, tryr, trz, trzl, trzr, vom(1:3), vpx, &
                                  vpy, vpz, x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable :: fxl(:), fxr(:), fxt(:), fyl(:), fyr(:), fyt(:), fzl(:), fzr(:), &
                                  fzt(:), ggx(:), ggy(:), ggz(:), oxt(:), oyt(:), ozt(:), q0t(:), &
                                  q1t(:), q2t(:), q3t(:), rgdoxt(:), rgdoyt(:), rgdozt(:), &
                                  rgdvxt(:), rgdvyt(:), rgdvzt(:), rgdxxt(:), rgdyyt(:), &
                                  rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), yyt(:), zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms), oyt(1:config%mxatms), ozt(1:config%mxatms), Stat=fail(6))
    End If
    Allocate (ggx(1:rigid%max_list * rigid%max_rigid), &
              ggy(1:rigid%max_list * rigid%max_rigid), &
              ggz(1:rigid%max_list * rigid%max_rigid), Stat=fail(7))
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(8))
    Allocate (vxt(1:config%mxatms), vyt(1:config%mxatms), vzt(1:config%mxatms), Stat=fail(9))
    Allocate (fxt(1:config%mxatms), fyt(1:config%mxatms), fzt(1:config%mxatms), Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
              q1t(1:rigid%max_rigid), &
              q2t(1:rigid%max_rigid), &
              q3t(1:rigid%max_rigid), Stat=fail(11))
    Allocate (rgdxxt(1:rigid%max_rigid), &
              rgdyyt(1:rigid%max_rigid), &
              rgdzzt(1:rigid%max_rigid), Stat=fail(12))
    Allocate (rgdvxt(1:rigid%max_rigid), &
              rgdvyt(1:rigid%max_rigid), &
              rgdvzt(1:rigid%max_rigid), Stat=fail(13))
    Allocate (rgdoxt(1:rigid%max_rigid), &
              rgdoyt(1:rigid%max_rigid), &
              rgdozt(1:rigid%max_rigid), Stat=fail(14))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l1 allocation failure'
      Call error(0, message)
    End If

    If (thermo%newjob_1) Then
      thermo%newjob_1 = .false.

      ! set number of constraint+pmf shake iterations

      If (cons%megcon > 0 .or. pmf%megpmf > 0) thermo%mxkit = 1
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake

      ! thermo%unsafe positioning due to possibly locally shared RBs

      thermo%unsafe = (Any(domain%map == comm%idnode))
    End If

    ! set matms

    matms = config%nlast
    If (comm%mxnode == 1) matms = config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
      End If
    End If

    ! Get the RB particles vectors wrt the RB's COM

    krgd = 0
    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then
        Do jrgd = 1, lrgd
          krgd = krgd + 1

          i = rigid%index_local(jrgd, irgd) ! local index of particle/site

          ! COM distances

          ggx(krgd) = config%parts(i)%xxx - rigid%xxx(irgd)
          ggy(krgd) = config%parts(i)%yyy - rigid%yyy(irgd)
          ggz(krgd) = config%parts(i)%zzz - rigid%zzz(irgd)
        End Do
      End If
    End Do

    ! minimum image convention for bond vectors

    Call images(config%imcon, config%cell, krgd, ggx, ggy, ggz)

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    rstep = 1.0_wp / tstep

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

      ! store initial values

      Do i = 1, matms
        xxt(i) = config%parts(i)%xxx
        yyt(i) = config%parts(i)%yyy
        zzt(i) = config%parts(i)%zzz

        vxt(i) = config%vxx(i)
        vyt(i) = config%vyy(i)
        vzt(i) = config%vzz(i)

        fxt(i) = config%parts(i)%fxx
        fyt(i) = config%parts(i)%fyy
        fzt(i) = config%parts(i)%fzz
      End Do

      Do irgd = 1, rigid%n_types
        q0t(irgd) = rigid%q0(irgd)
        q1t(irgd) = rigid%q1(irgd)
        q2t(irgd) = rigid%q2(irgd)
        q3t(irgd) = rigid%q3(irgd)

        rgdxxt(irgd) = rigid%xxx(irgd)
        rgdyyt(irgd) = rigid%yyy(irgd)
        rgdzzt(irgd) = rigid%zzz(irgd)

        rgdvxt(irgd) = rigid%vxx(irgd)
        rgdvyt(irgd) = rigid%vyy(irgd)
        rgdvzt(irgd) = rigid%vzz(irgd)

        rgdoxt(irgd) = rigid%oxx(irgd)
        rgdoyt(irgd) = rigid%oyy(irgd)
        rgdozt(irgd) = rigid%ozz(irgd)
      End Do

      ! Set afresh Langevin random forces

      Allocate (fxr(1:config%mxatms), fyr(1:config%mxatms), fzr(1:config%mxatms), Stat=fail(1))
      Allocate (fxl(1:config%mxatms), fyl(1:config%mxatms), fzl(1:config%mxatms), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'nvt_l1 allocation failure+'
        Call error(0, message)
      End If

      Call langevin_forces(nstep, thermo%temp, tstep, thermo%chi, fxr, fyr, fzr, cshell, config, seed)
      Call langevin_forces(-nstep, thermo%temp, tstep, thermo%chi, fxl, fyl, fzl, cshell, config, seed)
      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, fxr, fyr, fzr, domain, comm)
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, fxl, fyl, fzl, domain, comm)
      End If

      100 Continue

      ! constraint virial and stress tensor

      If (cons%megcon > 0) Then
        stat%vircon = 0.0_wp
        stat%strcon = 0.0_wp
      End If

      ! PMF virial and stress tensor

      If (pmf%megpmf > 0) Then
        stat%virpmf = 0.0_wp
        stat%strpmf = 0.0_wp
      End If

      ! Create primitive scalers and adjust/increase timestep if need be
      ! when Cholesky factorisation is compromised

      t0 = Exp(-thermo%chi * tstep)
      t1 = (1.0_wp - t0) / (thermo%chi)
      t2 = (1.0_wp - t0**2) / (2 * thermo%chi)

      safe = .true.
      Do
        tmp = t1**2 / t2
        If (tstep - tmp >= zero_plus) Then
          If (.not. safe) Then
            Write (message, "('timestep increased due to impossibility of integration, new timestep is:',3x,1p,e16.8,/)") &
              tstep
            Call info(message, .true.)
          End If
          Exit
        Else
          safe = .false.
          tstep = tmp + 1.0e-10_wp
          t0 = Exp(-thermo%chi * tstep)
          t1 = (1.0_wp - t0) / (thermo%chi)
          t2 = (1.0_wp - t0**2) / (2 * thermo%chi)
        End If
      End Do

      ! Create complex scalers

      scr = (t1 - t2) / Sqrt(t2)
      scl = Sqrt(tstep - (t1**2) / t2)
      scv = Sqrt(t2)

      scr1 = (t1 - t2) / Sqrt(t2 * tstep) / thermo%chi
      scl1 = Sqrt(Max(1.0_wp - (t1**2) / (t2 * tstep), 0.0_wp)) / thermo%chi
      scv1 = Sqrt(t2 / tstep)

      ! update velocity and position of FPs

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%weight(i) > 1.0e-6_wp) Then

          ! Half-kick velocity

          tmp = hstep / config%weight(i)
          config%vxx(i) = vxt(i) + tmp * fxt(i)
          config%vyy(i) = vyt(i) + tmp * fyt(i)
          config%vzz(i) = vzt(i) + tmp * fzt(i)

          tmp = tstep / config%weight(i)

          ! Full time fluctuations on positions using half-kick velocity

          config%parts(i)%xxx = xxt(i) + config%vxx(i) * t1 + tmp * (fxr(i) * scr1 + fxl(i) * scl1)
          config%parts(i)%yyy = yyt(i) + config%vyy(i) * t1 + tmp * (fyr(i) * scr1 + fyl(i) * scl1)
          config%parts(i)%zzz = zzt(i) + config%vzz(i) * t1 + tmp * (fzr(i) * scr1 + fzl(i) * scl1)

          ! Full time fluctuations on half-kick velocity

          config%vxx(i) = config%vxx(i) * t0 + tmp * fxr(i) * scv1
          config%vyy(i) = config%vyy(i) * t0 + tmp * fyr(i) * scv1
          config%vzz(i) = config%vzz(i) * t0 + tmp * fzr(i) * scv1

        End If
      End Do

      ! SHAKE procedures

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep, oxt, oyt, ozt, &
                         lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
      End If

      ! update velocity and position of RBs

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then

          ! calculate COM force and torque

          fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
          tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp

          fmxr = 0.0_wp; fmyr = 0.0_wp; fmzr = 0.0_wp
          tqxr = 0.0_wp; tqyr = 0.0_wp; tqzr = 0.0_wp

          fmxl = 0.0_wp; fmyl = 0.0_wp; fmzl = 0.0_wp
          tqxl = 0.0_wp; tqyl = 0.0_wp; tqzl = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + fxt(i)
              fmy = fmy + fyt(i)
              fmz = fmz + fzt(i)

              fmxr = fmxr + fxr(i)
              fmyr = fmyr + fyr(i)
              fmzr = fmzr + fzr(i)

              fmxl = fmxl + fxl(i)
              fmyl = fmyl + fyl(i)
              fmzl = fmzl + fzl(i)
            End If

            tqx = tqx + ggy(krgd) * fzt(i) - ggz(krgd) * fyt(i)
            tqy = tqy + ggz(krgd) * fxt(i) - ggx(krgd) * fzt(i)
            tqz = tqz + ggx(krgd) * fyt(i) - ggy(krgd) * fxt(i)

            tqxr = tqxr + ggy(krgd) * fzr(i) - ggz(krgd) * fyr(i)
            tqyr = tqyr + ggz(krgd) * fxr(i) - ggx(krgd) * fzr(i)
            tqzr = tqzr + ggx(krgd) * fyr(i) - ggy(krgd) * fxr(i)

            tqxl = tqxl + ggy(krgd) * fzl(i) - ggz(krgd) * fyl(i)
            tqyl = tqyl + ggz(krgd) * fxl(i) - ggx(krgd) * fzl(i)
            tqzl = tqzl + ggx(krgd) * fyl(i) - ggy(krgd) * fxl(i)
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation

          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = xxt(i1) - xxt(i2)
            y(1) = yyt(i1) - yyt(i2)
            z(1) = zzt(i1) - zzt(i2)

            Call images(config%imcon, config%cell, 1, x, y, z)

            mxdr = 1.0_wp / (x(1)**2 + y(1)**2 + z(1)**2)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) * mxdr
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp

            tmp = (x(1) * tqxr + y(1) * tqyr + z(1) * tqzr) * mxdr
            tqxr = x(1) * tmp
            tqyr = y(1) * tmp
            tqzr = z(1) * tmp

            tmp = (x(1) * tqxl + y(1) * tqyl + z(1) * tqzl) * mxdr
            tqxl = x(1) * tmp
            tqyl = y(1) * tmp
            tqzl = z(1) * tmp
          End If

          ! current rotation matrix

          Call getrotmat(q0t(irgd), q1t(irgd), q2t(irgd), q3t(irgd), rot)

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          trxr = tqxr * rot(1) + tqyr * rot(4) + tqzr * rot(7)
          tryr = tqxr * rot(2) + tqyr * rot(5) + tqzr * rot(8)
          trzr = tqxr * rot(3) + tqyr * rot(6) + tqzr * rot(9)

          trxl = tqxl * rot(1) + tqyl * rot(4) + tqzl * rot(7)
          tryl = tqxl * rot(2) + tqyl * rot(5) + tqzl * rot(8)
          trzl = tqxl * rot(3) + tqyl * rot(6) + tqzl * rot(9)

          ! calculate quaternion torques

          qt0 = 2.0_wp * (-q1t(irgd) * trx - q2t(irgd) * try - q3t(irgd) * trz)
          qt1 = 2.0_wp * (q0t(irgd) * trx - q3t(irgd) * try + q2t(irgd) * trz)
          qt2 = 2.0_wp * (q3t(irgd) * trx + q0t(irgd) * try - q1t(irgd) * trz)
          qt3 = 2.0_wp * (-q2t(irgd) * trx + q1t(irgd) * try + q0t(irgd) * trz)

          qt0r = 2.0_wp * (-q1t(irgd) * trxr - q2t(irgd) * tryr - q3t(irgd) * trzr)
          qt1r = 2.0_wp * (q0t(irgd) * trxr - q3t(irgd) * tryr + q2t(irgd) * trzr)
          qt2r = 2.0_wp * (q3t(irgd) * trxr + q0t(irgd) * tryr - q1t(irgd) * trzr)
          qt3r = 2.0_wp * (-q2t(irgd) * trxr + q1t(irgd) * tryr + q0t(irgd) * trzr)

          qt0l = 2.0_wp * (-q1t(irgd) * trxl - q2t(irgd) * tryl - q3t(irgd) * trzl)
          qt1l = 2.0_wp * (q0t(irgd) * trxl - q3t(irgd) * tryl + q2t(irgd) * trzl)
          qt2l = 2.0_wp * (q3t(irgd) * trxl + q0t(irgd) * tryl - q1t(irgd) * trzl)
          qt3l = 2.0_wp * (-q2t(irgd) * trxl + q1t(irgd) * tryl + q0t(irgd) * trzl)

          ! recover quaternion momenta at start of time step

          opx = rgdoxt(irgd) * rigid%rix(1, rgdtyp)
          opy = rgdoyt(irgd) * rigid%riy(1, rgdtyp)
          opz = rgdozt(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-q1t(irgd) * opx - q2t(irgd) * opy - q3t(irgd) * opz)
          p1 = 2.0_wp * (q0t(irgd) * opx - q3t(irgd) * opy + q2t(irgd) * opz)
          p2 = 2.0_wp * (q3t(irgd) * opx + q0t(irgd) * opy - q1t(irgd) * opz)
          p3 = 2.0_wp * (-q2t(irgd) * opx + q1t(irgd) * opy + q0t(irgd) * opz)

          ! update quaternion momenta to half-kick

          p0 = p0 + hstep * qt0
          p1 = p1 + hstep * qt1
          p2 = p2 + hstep * qt2
          p3 = p3 + hstep * qt3

          ! Full time fluctuations on quaternions

          tmp = t1 / tstep
          p00 = p0 * tmp + qt0r * scr1 + qt0l * scl1
          p11 = p1 * tmp + qt1r * scr1 + qt1l * scl1
          p22 = p2 * tmp + qt2r * scr1 + qt2l * scl1
          p33 = p3 * tmp + qt3r * scr1 + qt3l * scl1

          ! rotate RB quaternions - update q to full timestep & amend p
          ! and get new rotation matrix

          Call no_squish &
            (tstep, rigid%rix(2, rgdtyp), rigid%riy(2, rgdtyp), rigid%riz(2, rgdtyp), &
             rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), p00, p11, p22, p33)
          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! Full time fluctuations on quaternions - corrected

          p0 = (p00 - qt0r * scr1 - qt0l * scl1) / tmp
          p1 = (p11 - qt1r * scr1 - qt1l * scl1) / tmp
          p2 = (p22 - qt2r * scr1 - qt2l * scl1) / tmp
          p3 = (p33 - qt3r * scr1 - qt3l * scl1) / tmp

          ! Full time fluctuations on half-kick momenta

          p0 = p0 * t0 + tstep * qt0r * scv1
          p1 = p1 * t0 + tstep * qt1r * scv1
          p2 = p2 * t0 + tstep * qt2r * scv1
          p3 = p3 * t0 + tstep * qt3r * scv1

          ! update RB angular velocity to half step

          opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
          opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
          opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

          rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
          rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
          rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

          ! update RB COM velocity to half-kick

          tmp = hstep / rigid%weight(0, rgdtyp)
          rigid%vxx(irgd) = rgdvxt(irgd) + tmp * fmx
          rigid%vyy(irgd) = rgdvyt(irgd) + tmp * fmy
          rigid%vzz(irgd) = rgdvzt(irgd) + tmp * fmz

          tmp = tstep / rigid%weight(0, rgdtyp)

          ! Full time fluctuations on positions using half-kick velocity

          rigid%xxx(irgd) = rgdxxt(irgd) + rigid%vxx(irgd) * t1 + tmp * (fmxr * scr1 + fmxl * scl1)
          rigid%yyy(irgd) = rgdyyt(irgd) + rigid%vyy(irgd) * t1 + tmp * (fmyr * scr1 + fmyl * scl1)
          rigid%zzz(irgd) = rgdzzt(irgd) + rigid%vzz(irgd) * t1 + tmp * (fmzr * scr1 + fmzl * scl1)

          ! Full time fluctuations on half-kick velocity

          rigid%vxx(irgd) = rigid%vxx(irgd) * t0 + tmp * fmxr * scv1
          rigid%vyy(irgd) = rigid%vyy(irgd) * t0 + tmp * fmyr * scv1
          rigid%vzz(irgd) = rigid%vzz(irgd) * t0 + tmp * fmzr * scv1

          ! update RB members positions and halfstep velocities

          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              If (rigid%frozen(jrgd, rgdtyp) == 0) Then
                x(1) = rigid%x(jrgd, rgdtyp)
                y(1) = rigid%y(jrgd, rgdtyp)
                z(1) = rigid%z(jrgd, rgdtyp)

                ! new atomic positions

                config%parts(i)%xxx = rot(1) * x(1) + rot(2) * y(1) + rot(3) * z(1) + rigid%xxx(irgd)
                config%parts(i)%yyy = rot(4) * x(1) + rot(5) * y(1) + rot(6) * z(1) + rigid%yyy(irgd)
                config%parts(i)%zzz = rot(7) * x(1) + rot(8) * y(1) + rot(9) * z(1) + rigid%zzz(irgd)

                ! new atomic velocities in body frame

                vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                ! DD bound positions

                If (thermo%unsafe) Then
                  x(1) = config%parts(i)%xxx - xxt(i)
                  y(1) = config%parts(i)%yyy - yyt(i)
                  z(1) = config%parts(i)%zzz - zzt(i)
                  Call images(config%imcon, config%cell, 1, x, y, z)
                  config%parts(i)%xxx = x(1) + xxt(i)
                  config%parts(i)%yyy = y(1) + yyt(i)
                  config%parts(i)%zzz = z(1) + zzt(i)
                End If

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              End If
            End If
          End Do

        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, tmp, comm)) Then
          Call info(message, .true.)

          ! scale Langevin random forces

          Do i = 1, matms
            fxr(i) = fxr(i) * tmp
            fyr(i) = fyr(i) * tmp
            fzr(i) = fzr(i) * tmp

            fxl(i) = fxl(i) * tmp
            fyl(i) = fyl(i) * tmp
            fzl(i) = fzl(i) * tmp
          End Do

          ! restore initial conditions

          Do irgd = 1, rigid%n_types
            rigid%q0(irgd) = q0t(irgd)
            rigid%q1(irgd) = q1t(irgd)
            rigid%q2(irgd) = q2t(irgd)
            rigid%q3(irgd) = q3t(irgd)
          End Do

          ! restart vv1

          Go To 100
        End If
      End If

      Deallocate (fxr, fyr, fzr, Stat=fail(1))
      Deallocate (fxl, fyl, fzl, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'nvt_l1 deallocation failure+'
        Call error(0, message)
      End If

      ! second stage of velocity verlet algorithm
    Else

      ! update velocity (another half-kick) of FPs
      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
          config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
          config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! Get RB COM stress and virial
      Call rigid_bodies_stress(strcom, ggx, ggy, ggz, rigid, config, comm)
      vircom = -(strcom(1) + strcom(5) + strcom(9))

      ! update velocity of RBs
      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs
        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then ! Not that it matters

          ! calculate COM force and torque
          fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
          tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force
            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + config%parts(i)%fxx
              fmy = fmy + config%parts(i)%fyy
              fmz = fmz + config%parts(i)%fzz
            End If

            tqx = tqx + ggy(krgd) * config%parts(i)%fzz - ggz(krgd) * config%parts(i)%fyy
            tqy = tqy + ggz(krgd) * config%parts(i)%fxx - ggx(krgd) * config%parts(i)%fzz
            tqz = tqz + ggx(krgd) * config%parts(i)%fyy - ggy(krgd) * config%parts(i)%fxx
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation
          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = config%parts(i1)%xxx - config%parts(i2)%xxx
            y(1) = config%parts(i1)%yyy - config%parts(i2)%yyy
            z(1) = config%parts(i1)%zzz - config%parts(i2)%zzz

            Call images(config%imcon, config%cell, 1, x, y, z)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) / (x(1)**2 + y(1)**2 + z(1)**2)
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp
          End If

          ! current rotation matrix
          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! calculate torque in principal frame
          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          ! calculate quaternion torques
          qt0 = 2.0_wp * (-rigid%q1(irgd) * trx - rigid%q2(irgd) * try - rigid%q3(irgd) * trz)
          qt1 = 2.0_wp * (rigid%q0(irgd) * trx - rigid%q3(irgd) * try + rigid%q2(irgd) * trz)
          qt2 = 2.0_wp * (rigid%q3(irgd) * trx + rigid%q0(irgd) * try - rigid%q1(irgd) * trz)
          qt3 = 2.0_wp * (-rigid%q2(irgd) * trx + rigid%q1(irgd) * try + rigid%q0(irgd) * trz)

          ! recover quaternion momenta at half time step
          opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
          opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
          opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-rigid%q1(irgd) * opx - rigid%q2(irgd) * opy - rigid%q3(irgd) * opz)
          p1 = 2.0_wp * (rigid%q0(irgd) * opx - rigid%q3(irgd) * opy + rigid%q2(irgd) * opz)
          p2 = 2.0_wp * (rigid%q3(irgd) * opx + rigid%q0(irgd) * opy - rigid%q1(irgd) * opz)
          p3 = 2.0_wp * (-rigid%q2(irgd) * opx + rigid%q1(irgd) * opy + rigid%q0(irgd) * opz)

          ! update quaternion momenta to full step

          p0 = p0 + hstep * qt0
          p1 = p1 + hstep * qt1
          p2 = p2 + hstep * qt2
          p3 = p3 + hstep * qt3

          ! update RB angular & COM velocities to full step

          opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
          opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
          opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

          rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
          rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
          rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

          tmp = hstep / rigid%weight(0, rgdtyp)
          rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * fmx
          rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * fmy
          rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * fmz

          ! update RB members velocities

          Do jrgd = 1, lrgd
            If (rigid%frozen(jrgd, rgdtyp) == 0) Then
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site

              If (i <= config%natms) Then
                x(1) = rigid%x(jrgd, rgdtyp)
                y(1) = rigid%y(jrgd, rgdtyp)
                z(1) = rigid%z(jrgd, rgdtyp)

                ! new atomic velocities in body frame

                vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              End If
            End If
          End Do

        End If
      End Do

      ! remove system centre of mass velocity

      Call getvom(vom, rigid, config, comm)

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        If (rigid%frozen(0, rgdtyp) == 0) Then
          rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
          rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
          rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

          lrgd = rigid%list(-1, irgd)
          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              config%vxx(i) = config%vxx(i) - vom(1)
              config%vyy(i) = config%vyy(i) - vom(2)
              config%vzz(i) = config%vzz(i) - vom(3)
            End If
          End Do
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstresf(strknf, config, comm)
      Call kinstrest(rigid, strknt, comm)

      strkin = strknf + strknt
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

      ! update rotational energy

      engrot = getknr(rigid, comm)

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%deallocate_work()
      Deallocate (oxt, oyt, ozt, Stat=fail(6))
    End If
    Deallocate (ggx, ggy, ggz, Stat=fail(7))
    Deallocate (xxt, yyt, zzt, Stat=fail(8))
    Deallocate (vxt, vyt, vzt, Stat=fail(9))
    Deallocate (fxt, fyt, fzt, Stat=fail(10))
    Deallocate (q0t, q1t, q2t, q3t, Stat=fail(11))
    Deallocate (rgdxxt, rgdyyt, rgdzzt, Stat=fail(12))
    Deallocate (rgdvxt, rgdvyt, rgdvzt, Stat=fail(13))
    Deallocate (rgdoxt, rgdoyt, rgdozt, Stat=fail(14))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l1 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nvt_l1_vv

  Subroutine nvt_l2_vv(stage, lvar, mndis, mxdis, mxstp, tstep, nstep, strkin, engke, &
                       ttm, cshell, cons, pmf, stat, thermo, domain, tmr, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with inhomogeneous
    ! (two-temperature) Langevin thermostat, based on Langevin impulse
    ! (LI) integration
    !
    ! Ref: Jesus A. Izaguirre, `Langevin Stabilisation of Multiscale
    !      Mollified Dynamics', editors: A. Brandt, K. Binder, J. Bernholc,
    !      Multiscale Computational Methods in Chemistry and Physics,
    !      January 2001, vol. 117 of NATO Science Series: Series III -
    !      Computer and System Sciences, pp. 34-47, IOS Press, Amsterdam
    !
    ! (brownian dynamics is not symplectic due to the random forces)
    !
    ! copyright - daresbury laboratory
    ! authors   - i.t.todorov & m.a.seaton june 2017
    ! contrib   - m.a.seaton february 2020
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), engke
    Type(ttm_type),           Intent(InOut) :: ttm
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:9), i, ia, ijk, ja, ka
    Logical                    :: lrand, lvel, safe
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: chi, eltempmax, eltempmin, hstep, rstep, scl1, scl1a, scl1b, scr1, scr1a, &
                                  scr1b, scv1, scv1a, scv1b, t0, t0a, t0b, t1, t1a, t1b, t2, t2a, &
                                  t2b, tmp, velsq, vom(1:3)
    Real(Kind=wp), Allocatable :: fxl(:), fxr(:), fxt(:), fyl(:), fyr(:), fyt(:), fzl(:), fzr(:), &
                                  fzt(:), oxt(:), oyt(:), ozt(:), vxt(:), vyt(:), vzt(:), xxt(:), &
                                  yyt(:), zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms), oyt(1:config%mxatms), ozt(1:config%mxatms), Stat=fail(6))
    End If
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(7))
    Allocate (vxt(1:config%mxatms), vyt(1:config%mxatms), vzt(1:config%mxatms), Stat=fail(8))
    Allocate (fxt(1:config%mxatms), fyt(1:config%mxatms), fzt(1:config%mxatms), Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l2 allocation failure'
      Call error(0, message)
    End If

    If (thermo%newjob_2) Then
      thermo%newjob_2 = .false.

      ! set number of constraint+pmf shake iterations
      If (cons%megcon > 0 .or. pmf%megpmf > 0) thermo%mxkit = 1
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms
      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms
      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
      End If
    End If

    ! timestep derivatives
    hstep = 0.5_wp * tstep
    rstep = 1.0_wp / tstep

    ! Rescale chi to match average electronic temperature if
    ! using homogeneous electron-phonon coupling
    If (ttm%l_ttm .and. ttm%gvar == 1) Then
      Call calcchies(thermo%chi_ep, ttm, comm)
    End If

    ! check whether or not Langevin forces are needed: if electron-phonon
    ! friction coefficient is/will be greater than zero and coupling is
    ! switched on after time offset
    lrand = ((thermo%chi_ep > zero_plus .or. ttm%gvar == 2) .and. ttm%l_epcp)

    ! first pass of velocity verlet algorithm
    If (stage == VV_FIRST_STAGE) Then

      ! store initial values
      Do i = 1, config%natms
        xxt(i) = config%parts(i)%xxx
        yyt(i) = config%parts(i)%yyy
        zzt(i) = config%parts(i)%zzz

        vxt(i) = config%vxx(i)
        vyt(i) = config%vyy(i)
        vzt(i) = config%vzz(i)

        fxt(i) = config%parts(i)%fxx
        fyt(i) = config%parts(i)%fyy
        fzt(i) = config%parts(i)%fzz
      End Do

      ! Set afresh Langevin random forces
      Allocate (fxr(1:config%mxatms), fyr(1:config%mxatms), fzr(1:config%mxatms), Stat=fail(1))
      Allocate (fxl(1:config%mxatms), fyl(1:config%mxatms), fzl(1:config%mxatms), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'nvt_l2 allocation failure+'
        Call error(0, message)
      End If

      If (lrand) Then
        Call langevin_forces(nstep, thermo%temp, tstep, thermo%chi_ep, fxr, fyr, fzr, cshell, config, seed, ttm)
        Call langevin_forces(-nstep, thermo%temp, tstep, thermo%chi_ep, fxl, fyl, fzl, cshell, config, seed, ttm)
      Else
        fxr = 0.0_wp; fyr = 0.0_wp; fzr = 0.0_wp
        fxl = 0.0_wp; fyl = 0.0_wp; fzl = 0.0_wp
      End If

      100 Continue

      ! constraint virial and stress tensor

      If (cons%megcon > 0) Then
        stat%vircon = 0.0_wp
        stat%strcon = 0.0_wp
      End If

      ! PMF virial and stress tensor

      If (pmf%megpmf > 0) Then
        stat%virpmf = 0.0_wp
        stat%strpmf = 0.0_wp
      End If

      ! Create primitive scalers and adjust/increase timestep if need be
      ! when Cholesky factorisation is compromised (using minimum possible
      ! value of friction factor to give largest timestep)

      Select Case (ttm%gvar)
      Case (0, 1)
        chi = Min(thermo%chi_ep, thermo%chi_ep + thermo%chi_es)
      Case (2)
        Call eltemp_max(eltempmax, ttm, comm)
        Call eltemp_min(eltempmin, ttm, comm)
        chi = Gep(eltempmin, ttm)
        chi = Min(chi, Gep(eltempmax, ttm))
        chi = Min(chi, chi + thermo%chi_es)
      End Select
      t0 = Exp(-chi * tstep)
      t1 = (1.0_wp - t0) / (chi)
      t2 = (1.0_wp - t0**2) / (2 * chi)

      safe = .true.
      Do
        tmp = t1**2 / t2
        If (tstep - tmp >= zero_plus) Then
          If (.not. safe) Then
            Write (message, "('timestep increased due to impossibility of integration, new timestep is:',3x,1p,e16.8)") &
              tstep
            Call info(message, .true.)
          End If
          Exit
        Else
          safe = .false.
          tstep = tmp + 1.0e-10_wp
          t0 = Exp(-chi * tstep)
          t1 = (1.0_wp - t0) / (chi)
          t2 = (1.0_wp - t0**2) / (2 * chi)
        End If
      End Do

      ! Create complex scalers: only for constant and homogeneous
      ! electron-phonon coupling

      Select Case (ttm%gvar)
      Case (0, 1)
        chi = Merge(thermo%chi_ep, 0.0_wp, ttm%l_epcp) + thermo%chi_es
        t0a = Exp(-chi * tstep)
        If (chi > zero_plus) Then
          t1a = (1.0_wp - t0a) / (chi)
          t2a = (1.0_wp - t0a**2) / (2 * chi)
          scr1a = (t1a - t2a) / Sqrt(t2a * tstep) / chi
          scl1a = Sqrt(1.0_wp - (t1a**2) / (t2a * tstep)) / chi
          scv1a = Sqrt(t2a / tstep)
        Else
          t1a = tstep
          t2a = 0.0_wp
          scr1a = 0.0_wp
          scl1a = 0.0_wp
          scv1a = 0.0_wp
        End If
        If (lrand) Then
          t0b = Exp(-thermo%chi_ep * tstep)
          t1b = (1.0_wp - t0b) / thermo%chi_ep
          t2b = (1.0_wp - t0b**2) / (2 * thermo%chi_ep)
          scr1b = (t1b - t2b) / Sqrt(t2b * tstep) / thermo%chi_ep
          scl1b = Sqrt(1.0_wp - (t1b**2) / (t2b * tstep)) / thermo%chi_ep
          scv1b = Sqrt(t2b / tstep)
        Else
          t0b = 0.0_wp
          t1b = tstep
          t2b = 0.0_wp
          scr1b = 0.0_wp
          scl1b = 0.0_wp
          scv1b = 0.0_wp
        End If
      End Select

      ! update velocity and position

      If (ttm%l_ttm) Then

        If (ttm%oneway) Then
          ! one-way electron-phonon coupling
          Do i = 1, config%natms
            velsq = config%vxx(i) * config%vxx(i) + config%vyy(i) * config%vyy(i) + config%vzz(i) * config%vzz(i)
            lvel = (velsq > thermo%vel_es2 .and. thermo%chi_es > zero_plus)
            If (config%weight(i) > 1.0e-6_wp) Then
              ! check for active config%cell and electronic temperature is
              ! higher than ionic tmeperature: if not, switch off thermostat
              ia = Floor(xxt(i)*ttm%grcell(1)+yyt(i)*ttm%grcell(4)+zzt(i)*ttm%grcell(7)+ttm%zerocell(1)) + 1
              ja = Floor(xxt(i)*ttm%grcell(2)+yyt(i)*ttm%grcell(5)+zzt(i)*ttm%grcell(8)+ttm%zerocell(2)) + 1
              ka = Floor(xxt(i)*ttm%grcell(3)+yyt(i)*ttm%grcell(6)+zzt(i)*ttm%grcell(9)+ttm%zerocell(3)) + 1
              ijk = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
              If (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus .and. ttm%eltemp(ijk, 0, 0, 0) > ttm%tempion(ijk)) Then
                Select Case (ttm%gvar)
                Case (0, 1)
                  t0 = Merge(t0a, t0b, lvel)
                  t1 = Merge(t1a, t1b, lvel)
                  scr1 = Merge(scr1a, scr1b, lvel)
                  scl1 = Merge(scl1a, scl1b, lvel)
                  scv1 = Merge(scv1a, scv1b, lvel)
                Case (2)
                  chi = Merge(Gep(ttm%eltemp(ijk, 0, 0, 0), ttm), 0.0_wp, ttm%l_epcp) + Merge(thermo%chi_es, 0.0_wp, lvel)
                  If (ttm%l_epcp) Then
                    t0 = Exp(-tstep * chi)
                    t1 = (1.0_wp - t0) / chi
                    t2 = (1.0_wp - t0**2) / (2 * chi)
                    scr1 = (t1 - t2) / Sqrt(t2 * tstep) / chi
                    scl1 = Sqrt(1.0_wp - (t1**2) / (t2 * tstep)) / chi
                    scv1 = Sqrt(t2 / tstep)
                  Else
                    t0 = 1.0_wp
                    t1 = tstep
                    scr1 = 0.0_wp
                    scl1 = 0.0_wp
                    scv1 = 0.0_wp
                  End If
                End Select
              Else
                t0 = 1.0_wp
                t1 = tstep
                scr1 = 0.0_wp
                scl1 = 0.0_wp
                scv1 = 0.0_wp
              End If

              ! Half-kick velocity

              tmp = hstep / config%weight(i)
              config%vxx(i) = vxt(i) + tmp * fxt(i)
              config%vyy(i) = vyt(i) + tmp * fyt(i)
              config%vzz(i) = vzt(i) + tmp * fzt(i)

              tmp = tstep / config%weight(i)

              ! Full time fluctuations on positions using half-kick velocity

              config%parts(i)%xxx = xxt(i) + config%vxx(i) * t1 + tmp * (fxr(i) * scr1 + fxl(i) * scl1)
              config%parts(i)%yyy = yyt(i) + config%vyy(i) * t1 + tmp * (fyr(i) * scr1 + fyl(i) * scl1)
              config%parts(i)%zzz = zzt(i) + config%vzz(i) * t1 + tmp * (fzr(i) * scr1 + fzl(i) * scl1)

              ! Full time fluctuations on half-kick perculiar (thermal only) velocity

              config%vxx(i) = (config%vxx(i) - ttm%ttmvom(ijk, 1)) * t0 + tmp * fxr(i) * scv1 + ttm%ttmvom(ijk, 1)
              config%vyy(i) = (config%vyy(i) - ttm%ttmvom(ijk, 2)) * t0 + tmp * fyr(i) * scv1 + ttm%ttmvom(ijk, 2)
              config%vzz(i) = (config%vzz(i) - ttm%ttmvom(ijk, 3)) * t0 + tmp * fzr(i) * scv1 + ttm%ttmvom(ijk, 3)

            End If
          End Do

        Else

          Do i = 1, config%natms
            velsq = config%vxx(i) * config%vxx(i) + config%vyy(i) * config%vyy(i) + config%vzz(i) * config%vzz(i)
            lvel = (velsq > thermo%vel_es2 .and. thermo%chi_es > zero_plus)
            If (config%weight(i) > 1.0e-6_wp) Then
              ! check for active config%cell: if not, switch off thermostat
              ia = Floor(xxt(i)*ttm%grcell(1)+yyt(i)*ttm%grcell(4)+zzt(i)*ttm%grcell(7)+ttm%zerocell(1)) + 1
              ja = Floor(xxt(i)*ttm%grcell(2)+yyt(i)*ttm%grcell(5)+zzt(i)*ttm%grcell(8)+ttm%zerocell(2)) + 1
              ka = Floor(xxt(i)*ttm%grcell(3)+yyt(i)*ttm%grcell(6)+zzt(i)*ttm%grcell(9)+ttm%zerocell(3)) + 1
              ijk = 1 + ia + (ttm%ntcell(1) + 2) * (ja + (ttm%ntcell(2) + 2) * ka)
              If (ttm%act_ele_cell(ijk, 0, 0, 0) > zero_plus) Then
                Select Case (ttm%gvar)
                Case (0, 1)
                  t0 = Merge(t0a, t0b, lvel)
                  t1 = Merge(t1a, t1b, lvel)
                  scr1 = Merge(scr1a, scr1b, lvel)
                  scl1 = Merge(scl1a, scl1b, lvel)
                  scv1 = Merge(scv1a, scv1b, lvel)
                Case (2)
                  chi = Merge(Gep(ttm%eltemp(ijk, 0, 0, 0), ttm), 0.0_wp, ttm%l_epcp) + Merge(thermo%chi_es, 0.0_wp, lvel)
                  If (ttm%l_epcp) Then
                    t0 = Exp(-tstep * chi)
                    t1 = (1.0_wp - t0) / chi
                    t2 = (1.0_wp - t0**2) / (2 * chi)
                    scr1 = (t1 - t2) / Sqrt(t2 * tstep) / chi
                    scl1 = Sqrt(1.0_wp - (t1**2) / (t2 * tstep)) / chi
                    scv1 = Sqrt(t2 / tstep)
                  Else
                    t0 = 1.0_wp
                    t1 = tstep
                    t2 = 0.0_wp
                    scr1 = 0.0_wp
                    scl1 = 0.0_wp
                    scv1 = 0.0_wp
                  End If
                End Select
              Else
                t0 = 1.0_wp
                t1 = tstep
                t2 = 0.0_wp
                scr1 = 0.0_wp
                scl1 = 0.0_wp
                scv1 = 0.0_wp
              End If

              ! Half-kick velocity

              tmp = hstep / config%weight(i)
              config%vxx(i) = vxt(i) + tmp * fxt(i)
              config%vyy(i) = vyt(i) + tmp * fyt(i)
              config%vzz(i) = vzt(i) + tmp * fzt(i)

              tmp = tstep / config%weight(i)

              ! Full time fluctuations on positions using half-kick velocity

              config%parts(i)%xxx = xxt(i) + config%vxx(i) * t1 + tmp * (fxr(i) * scr1 + fxl(i) * scl1)
              config%parts(i)%yyy = yyt(i) + config%vyy(i) * t1 + tmp * (fyr(i) * scr1 + fyl(i) * scl1)
              config%parts(i)%zzz = zzt(i) + config%vzz(i) * t1 + tmp * (fzr(i) * scr1 + fzl(i) * scl1)

              ! Full time fluctuations on half-kick peculiar (thermal only) velocity

              config%vxx(i) = (config%vxx(i) - ttm%ttmvom(ijk, 1)) * t0 + tmp * fxr(i) * scv1 + ttm%ttmvom(ijk, 1)
              config%vyy(i) = (config%vyy(i) - ttm%ttmvom(ijk, 2)) * t0 + tmp * fyr(i) * scv1 + ttm%ttmvom(ijk, 2)
              config%vzz(i) = (config%vzz(i) - ttm%ttmvom(ijk, 3)) * t0 + tmp * fzr(i) * scv1 + ttm%ttmvom(ijk, 3)

            End If
          End Do
        End If

      Else

        ! no ttm option: just inhomogeneous Langevin thermostat

        Do i = 1, config%natms
          velsq = config%vxx(i) * config%vxx(i) + config%vyy(i) * config%vyy(i) + config%vzz(i) * config%vzz(i)
          lvel = (velsq > thermo%vel_es2 .and. thermo%chi_es > zero_plus)
          If (config%weight(i) > 1.0e-6_wp) Then

            ! Half-kick velocity

            tmp = hstep / config%weight(i)
            config%vxx(i) = vxt(i) + tmp * fxt(i)
            config%vyy(i) = vyt(i) + tmp * fyt(i)
            config%vzz(i) = vzt(i) + tmp * fzt(i)

            t0 = Merge(t0a, t0b, lvel)
            t1 = Merge(t1a, t1b, lvel)
            scr1 = Merge(scr1a, scr1b, lvel)
            scl1 = Merge(scl1a, scl1b, lvel)
            scv1 = Merge(scv1a, scv1b, lvel)
            tmp = tstep / config%weight(i) * (1.0_wp + Merge(thermo%chi_es / thermo%chi_ep, 0.0_wp, lvel))

            ! Full time fluctuations on positions using half-kick velocity
            ! (time multipler adjusted to increase random forces when
            ! electron stopping is required due to increase in chi)

            config%parts(i)%xxx = xxt(i) + config%vxx(i) * t1 + tmp * (fxr(i) * scr1 + fxl(i) * scl1)
            config%parts(i)%yyy = yyt(i) + config%vyy(i) * t1 + tmp * (fyr(i) * scr1 + fyl(i) * scl1)
            config%parts(i)%zzz = zzt(i) + config%vzz(i) * t1 + tmp * (fzr(i) * scr1 + fzl(i) * scl1)

            ! Full time fluctuations on half-kick velocity

            config%vxx(i) = config%vxx(i) * t0 + tmp * fxr(i) * scv1
            config%vyy(i) = config%vyy(i) * t0 + tmp * fyr(i) * scv1
            config%vzz(i) = config%vzz(i) * t0 + tmp * fzr(i) * scv1

          End If
        End Do

      End If

      ! SHAKE procedures

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep, oxt, oyt, ozt, &
                         lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
      End If

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, tmp, comm)) Then
          Call info(message, .true.)

          ! scale Langevin random forces

          Do i = 1, config%natms
            fxr(i) = fxr(i) * tmp
            fyr(i) = fyr(i) * tmp
            fzr(i) = fzr(i) * tmp

            fxl(i) = fxl(i) * tmp
            fyl(i) = fyl(i) * tmp
            fzl(i) = fzl(i) * tmp
          End Do

          ! restart vv1

          Go To 100
        End If
      End If

      Deallocate (fxr, fyr, fzr, Stat=fail(1))
      Deallocate (fxl, fyl, fzl, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a,i0)') 'nvt_l2 deallocation failure+'
        Call error(0, message)
      End If

      ! second stage of velocity verlet algorithm

    Else

      ! update velocity (another half-kick)

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
          config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
          config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! remove system centre of mass velocity

      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstress(strkin, config, comm)
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%deallocate_work()
      Deallocate (oxt, oyt, ozt, Stat=fail(6))
    End If
    Deallocate (xxt, yyt, zzt, Stat=fail(7))
    Deallocate (vxt, vyt, vzt, Stat=fail(8))
    Deallocate (fxt, fyt, fzt, Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_l2 deallocation failure'
      Call error(0, message)
    End If
  End Subroutine nvt_l2_vv
End Module nvt_langevin
